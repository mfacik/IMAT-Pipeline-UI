import os
import sys
import shutil
import subprocess
from pathlib import Path

# Gurobi yerine açık kaynaklı GLPK çözücüsünü varsayılan yapıyoruz
os.environ["IMAT_SOLVER"] = "GLPK"

import pandas as pd
import numpy as np
import plotly.express as px
from shiny import App, render, reactive, ui

# ------------------------------------------------------------------------------
# 1. DOSYA YOLLARI VE DİZİN YAPISI
# ------------------------------------------------------------------------------
ROOT_DIR = Path(__file__).parent
RESOURCES_DIR = ROOT_DIR / "Resources"

# Demo Verileri Dizini
SHOWCASE_RESULTS = RESOURCES_DIR / "showcase" / "results"
SHOWCASE_FIGURES = RESOURCES_DIR / "showcase" / "figures"

# Canlı Pipeline Veri ve Çıktı Dizinleri
RAW_DATA_DIR = RESOURCES_DIR / "data" / "raw"
SYNTHETIC_COHORT_DIR = RESOURCES_DIR / "synthetic_single_cohort"
CANLI_RESULTS_DIR = RESOURCES_DIR / "results"

# ------------------------------------------------------------------------------
# 2. ARAYÜZ (UI)
# ------------------------------------------------------------------------------
app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h4("Kontrol Paneli"),
        
        ui.input_radio_buttons(
            "data_mode", 
            "Çalışma Modu:", 
            {
                "demo": "Demo Modu (Hazır MOCK Veri)", 
                "custom": "Kullanıcı Verisi Yükle"
            },
            selected="demo"
        ),
        
        ui.hr(),
        
        ui.input_radio_buttons(
            "selected_layer",
            "Katman Seçimi:",
            {
                "expression_only": "1. Layer: Expression Only",
                "splicing_informed": "2. Layer: Splicing Informed",
                "splicing_pathogenicity_informed": "3. Layer: Splicing + Pathogenicity"
            },
            selected="expression_only"
        ),
        
        ui.hr(),
        
        # KULLANICI VERİSİ YÜKLEME PANELSİ
        ui.panel_conditional(
            "input.data_mode == 'custom'",
            ui.h5("📁 Veri Yükleme & QC"),
            
            ui.input_file("file_expression", "1. Gen İfade Matrisi (.tsv):", accept=[".tsv", ".txt"]),
            ui.input_file("file_metadata", "2. Örnek Meta Verisi (.tsv):", accept=[".tsv", ".txt"]),
            ui.input_file("file_splicing", "3. Splicing Hedefleri (.tsv - İsteğe Bağlı):", accept=[".tsv", ".txt"]),
            ui.input_file("file_patho", "4. Patojenite Sıralaması (.tsv - İsteğe Bağlı):", accept=[".tsv", ".txt"]),
            
            ui.hr(),
            ui.h6("📋 Şablon İndir"),
            ui.download_button("download_template_expr", "1. Gen İfade Şablonu", class_="btn-sm btn-outline-secondary mb-1 w-100"),
            ui.download_button("download_template_meta", "2. Meta Veri Şablonu", class_="btn-sm btn-outline-secondary mb-1 w-100"),
            ui.download_button("download_template_splicing", "3. Splicing Şablonu", class_="btn-sm btn-outline-secondary mb-1 w-100"),
            ui.download_button("download_template_patho", "4. Patojenite Şablonu", class_="btn-sm btn-outline-secondary mb-1 w-100"),
            
            ui.hr()
        ),
        
        ui.panel_conditional(
            "input.selected_layer == 'splicing_pathogenicity_informed' || input.data_mode == 'custom'",
            ui.input_slider("patho_cutoff", "Patojenite Kesim Eşiği (Rank Threshold):", min=0.50, max=0.99, value=0.95, step=0.01)
        ),
        
        ui.input_slider("p_cutoff", "Max FDR/p-value Eşiği:", min=0.001, max=0.1, value=0.05, step=0.005),
        
        ui.hr(),
        
        ui.input_action_button("btn_start", "🚀 Analizi Başlat / Göster", class_="btn-primary w-100"),
        
        title="Parametreler"
    ),
    
    ui.output_ui("dynamic_tabs"),
    
    title="Splicing & Pathogenicity iMAT Dashboard",
    fillable=True
)

# ------------------------------------------------------------------------------
# 3. SUNUCU MANTIĞI (SERVER)
# ------------------------------------------------------------------------------
def server(input, output, session):

    log_messages = reactive.value("Terminal bekleniyor... Lütfen analizi başlatın.\n")
    progress_percent = reactive.value(0)
    progress_detail = reactive.value("Henüz işlem başlatılmadı.")
    
    analysis_started = reactive.value(False)
    pipeline_trigger = reactive.value(0)

    # --------------------------------------------------------------------------
    # DİNAMİK SEKMELERİN RENDER EDİLMESİ
    # --------------------------------------------------------------------------
    @render.ui
    def dynamic_tabs():
        about_panel = ui.nav_panel(
            "ℹ️ Genel Bakış (About)",
            ui.markdown("""
            ## Splicing & Pathogenicity-Informed iMAT Workbench
            Bu uygulama, gen ifade verileri (*Gene Expression*), alternatif eklenme (*Alternative Splicing*) hedefleri ve varyant patojenite skorlarını **Human-GEM 1.12.0** modeline entegre ederek dokusal/örneğe özgü reaksiyon aktivitelerini tahmin eder.
            """),
            value="tab_about"
        )

        if not analysis_started():
            return ui.navset_card_tab(about_panel, id="main_tabs")
        
        status_panel = ui.nav_panel(
            "⚙️ Pipeline Çalışma Durumu",
            ui.markdown("### Analiz Süreci ve İlerleme Takibi"),
            ui.card(
                ui.card_header("İşlem İlerleme Seviyesi"),
                ui.output_ui("progress_bar_ui")
            ),
            ui.card(
                ui.card_header("Canlı Terminal Çıktıları & QC Logları"),
                ui.tags.div(
                    ui.output_code("pipeline_logs"),
                    style="background-color: #1e1e1e; color: #ffffff; padding: 15px; border-radius: 5px; font-family: 'Courier New', monospace; max-height: 400px; overflow-y: auto;"
                )
            ),
            value="tab_status"
        )
        
        results_panel = ui.nav_panel(
            "📊 Analiz Sonuçları (Results)",
            ui.markdown("### Model Aktivite ve Büyüklük Dağılımları"),
            ui.layout_columns(
                ui.card(
                    ui.card_header("Reaksiyon Aktivite Dağılımı (Active Reactions)"),
                    ui.output_image("img_active_reactions")
                ),
                ui.card(
                    ui.card_header("Model Büyüklük Değişimi (Paired Model-Size Change)"),
                    ui.output_image("img_paired_size")
                )
            ),
            value="tab_results"
        )
        
        pathways_panel = ui.nav_panel(
            "🧬 Yolak Analizi (Pathways)",
            ui.markdown("### Anlamlı Değişen Metabolik Yolaklar (Reporter-like Analysis)"),
            ui.card(
                ui.card_header("Detaylı Yolak Tablosu"),
                ui.output_data_frame("table_pathways")
            ),
            value="tab_pathways"
        )

        return ui.navset_card_tab(
            about_panel,
            status_panel,
            results_panel,
            pathways_panel,
            id="main_tabs",
            selected="tab_results" if progress_percent() == 100 else "tab_status"
        )

    # --------------------------------------------------------------------------
    # DOSYALARI PİPELİNE'IN SIKINTISIZ OKUYACAĞI HER İKİ KLASÖRE DE YAZMA
    # --------------------------------------------------------------------------
    def save_uploaded_files_to_raw():
        RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)
        SYNTHETIC_COHORT_DIR.mkdir(parents=True, exist_ok=True)
        
        file_mappings = [
            (input.file_expression(), "expression.tsv"),
            (input.file_metadata(), "metadata.tsv"),
            (input.file_splicing(), "splicing_targets.tsv"),
            (input.file_patho(), "pathogenicity_ranks.tsv")
        ]
        
        for file_input, target_name in file_mappings:
            if file_input:
                src = file_input[0]["datapath"]
                shutil.copy(src, RAW_DATA_DIR / target_name)
                shutil.copy(src, SYNTHETIC_COHORT_DIR / target_name)

    def validate_custom_inputs():
        errors = []
        if input.file_expression() is None:
            errors.append("Eksik Dosya: Gen İfade Matrisi (.tsv) yüklenmedi.")
        if input.file_metadata() is None:
            errors.append("Eksik Dosya: Örnek Meta Verisi (.tsv) yüklenmedi.")
        else:
            try:
                meta_df = pd.read_csv(input.file_metadata()[0]["datapath"], sep="\t")
                if not {"sample_id", "group"}.issubset(meta_df.columns):
                    errors.append("Metadata Format Hatası: 'sample_id' ve 'group' sütunları bulunamadı.")
                else:
                    groups = set(meta_df["group"].unique())
                    if not groups.issubset({"ad", "control"}):
                        errors.append(f"Metadata İçerik Hatası: 'group' sütunu sadece 'ad' ve 'control' içermelidir.")
            except Exception as e:
                errors.append(f"Metadata Okuma Hatası: {str(e)}")

        return errors

    # --------------------------------------------------------------------------
    # ANALİZ BAŞLATMA EVENT HANDLER (ZİNCİRLEME SCRIPT ÇALIŞTIRMA)
    # --------------------------------------------------------------------------
    @reactive.effect
    @reactive.event(input.btn_start)
    def handle_analysis_start():
        mode = input.data_mode()
        
        if mode == "demo":
            analysis_started.set(True)
            progress_percent.set(30)
            progress_detail.set("Hazır Demo Verileri Yükleniyor...")
            log_messages.set(">> [DEMO MODU] Hazır MOCK veri seti seçildi.\n>> Veriler 'showcase/results' klasöründen yükleniyor...\n")
            
            progress_percent.set(100)
            progress_detail.set("Demo Verisi Başarıyla Yüklendi!")
            log_messages.set(log_messages() + ">> [BAŞARILI] İşlem tamamlandı. Sonuç sekmeleri aktifleştirildi.\n")
            ui.notification_show("Demo verisi başarıyla yüklendi!", type="default")
            
            pipeline_trigger.set(pipeline_trigger() + 1)
            
        elif mode == "custom":
            validation_errors = validate_custom_inputs()
            
            if validation_errors:
                analysis_started.set(True)
                progress_percent.set(0)
                progress_detail.set("Girdi Kontrolü Başarısız!")
                error_summary = "\n".join([f"❌ {err}" for err in validation_errors])
                log_messages.set(">> [GİRDİ KONTROLÜ BAŞARISIZ]\n" + error_summary + "\n")
                ui.notification_show("Lütfen eksik veya hatalı dosyaları düzeltin!", type="error")
                return

            analysis_started.set(True)
            save_uploaded_files_to_raw()
            log_messages.set(">> [GİRDİ KONTROLÜ BAŞARILI] Yüklenen dosyalar sisteme aktarıldı.\n")
            
            # --- YENİ EKLENEN KISIM: İLK KULLANIM İÇİN OTOMATİK MODEL İNDİRME ---
            model_file = RESOURCES_DIR / "model" / "Human-GEM-1.12.0.xml"
            if not model_file.exists():
                progress_percent.set(10)
                progress_detail.set("Sistem Hazırlığı: Referans Model İndiriliyor (Sadece İlk Kullanımda)...")
                log_messages.set(log_messages() + ">> [SİSTEM] Referans harita (Human-GEM-1.12.0) eksik. Otomatik olarak indiriliyor...\n")
                
                download_script = ROOT_DIR / "download_model.py"
                if not download_script.exists():
                    download_script = RESOURCES_DIR / "download_model.py"
                
                success = run_script(download_script)
                if not success:
                    progress_detail.set("Hata: Referans Model İndirilemedi!")
                    ui.notification_show("Model indirme başarısız oldu. İnternet bağlantınızı kontrol edin.", type="error")
                    return
                log_messages.set(log_messages() + ">> [BAŞARILI] Model başarıyla indirildi ve sisteme entegre edildi!\n")
            # ----------------------------------------------------------------------
            
            layer = input.selected_layer()
            
            # --- ZİNCİRLEME SCRIPT LİSTESİ (00_PREPARE ATLANDI!) ---
            # Kullanıcının kendi verisini ezmemesi için doğrudan 01'den başlıyoruz
            scripts_to_run = ["01_expression_only.py"]
            
            if layer == "splicing_informed":
                scripts_to_run.append("02_splicing_informed.py")
            elif layer == "splicing_pathogenicity_informed":
                scripts_to_run.append("02_splicing_informed.py")
                scripts_to_run.append("03_splicing_pathogenicity_informed.py")
                
            scripts_to_run.extend(["04_statistics_and_reporter_like.py", "05_plot_results.py"])
            
            total_steps = len(scripts_to_run)
            
            for idx, script_name in enumerate(scripts_to_run, start=1):
                pct = int((idx / total_steps) * 100)
                progress_percent.set(pct)
                progress_detail.set(f"Adım {idx}/{total_steps}: {script_name} çalıştırılıyor...")
                log_messages.set(log_messages() + f">> [{idx}/{total_steps}] {script_name} başlatılıyor...\n")
                
                success = run_script(RESOURCES_DIR / script_name)
                if not success:
                    progress_detail.set(f"Hata: {script_name} adımında işlem durduruldu!")
                    ui.notification_show(f"{script_name} çalışırken hata oluştu. Logları inceleyin.", type="error")
                    return
            
            progress_percent.set(100)
            progress_detail.set("Tüm Pipeline Adımları Başarıyla Tamamlandı!")
            log_messages.set(log_messages() + ">> [BAŞARILI] Gerçek verilerinizle yapılan analiz tamamlandı!\n")
            ui.notification_show("Kullanıcı verinizle analiz tamamlandı!", type="default")
            
            pipeline_trigger.set(pipeline_trigger() + 1)

    def run_script(script_path, extra_args=None):
        if extra_args is None:
            extra_args = []
        try:
            cmd = [sys.executable, str(script_path)] + extra_args
            # Eğer çalıştırılan script download_model ise ROOT_DIR'den çalıştır, yoksa RESOURCES_DIR'den.
            working_dir = str(ROOT_DIR) if "download_model" in script_path.name else str(RESOURCES_DIR)
            
            result = subprocess.run(
                cmd,
                cwd=working_dir,
                capture_output=True,
                text=True,
                check=True
            )
            log_messages.set(log_messages() + result.stdout + "\n")
            return True
        except subprocess.CalledProcessError as e:
            log_messages.set(log_messages() + f"\n[HATA] {script_path.name} çalışırken hata oluştu:\n" + e.stderr + "\n")
            return False

    # --------------------------------------------------------------------------
    # DİNAMİK VERİ VE GRAFİK OKUMA MANTIĞI (%100 GERÇEK VE CANLI)
    # --------------------------------------------------------------------------
    @reactive.calc
    def loaded_reporter_data():
        _ = pipeline_trigger()
        mode = input.data_mode()
        layer = input.selected_layer()
        
        if mode == "custom":
            file_path = CANLI_RESULTS_DIR / layer / "reporter_like" / "reporter_like_all_pathways.tsv"
            if not file_path.exists():
                file_path = CANLI_RESULTS_DIR / layer / "reporter_like_all_pathways.tsv"
            if not file_path.exists():
                return pd.DataFrame()
        else:
            file_path = SHOWCASE_RESULTS / layer / "reporter_like_all_pathways.tsv"
            if not file_path.exists():
                return pd.DataFrame()
                
        return pd.read_csv(file_path, sep="\t")

    @render.image
    def img_active_reactions():
        _ = pipeline_trigger()
        mode = input.data_mode()
        
        if mode == "custom":
            svg_path = CANLI_RESULTS_DIR / "figures" / "active_reactions.svg"
            if not svg_path.exists():
                return None
        else:
            svg_path = SHOWCASE_FIGURES / "active_reactions.svg"
            if not svg_path.exists():
                return None
                
        return {"src": str(svg_path), "width": "100%"}

    @render.image
    def img_paired_size():
        _ = pipeline_trigger()
        mode = input.data_mode()
        
        if mode == "custom":
            svg_path = CANLI_RESULTS_DIR / "figures" / "paired_model_size.svg"
            if not svg_path.exists():
                return None
        else:
            svg_path = SHOWCASE_FIGURES / "paired_model_size.svg"
            if not svg_path.exists():
                return None
                
        return {"src": str(svg_path), "width": "100%"}

    @render.data_frame
    def table_pathways():
        df = loaded_reporter_data()
        if df.empty:
            return render.DataGrid(pd.DataFrame())
        
        p_col = "reporter_p_corrected" if "reporter_p_corrected" in df.columns else "reporter_p_adjusted_bh"
        if p_col in df.columns:
            filtered = df[df[p_col] <= input.p_cutoff()]
        else:
            filtered = df
            
        cols = [c for c in ["subsystem", "reaction_count", "reporter_p_corrected", "reporter_p_adjusted_bh"] if c in filtered.columns]
        return render.DataGrid(filtered[cols], selection_mode="rows")

    # --------------------------------------------------------------------------
    # ŞABLON İNDİRME HANDLERLARI
    # --------------------------------------------------------------------------
    @render.download_button(filename="expression_template.tsv")
    def download_template_expr():
        yield "gene_id\tbiotype\tSample_001\tSample_002\nENSG00000141510\tprotein_coding\t15.4\t0.0\nENSG00000121410\tprotein_coding\t2.1\t45.8\n"

    @render.download_button(filename="metadata_template.tsv")
    def download_template_meta():
        yield "sample_id\tgroup\nSample_001\tad\nSample_002\tcontrol\n"

    @render.download_button(filename="splicing_targets_template.tsv")
    def download_template_splicing():
        yield "sample_id\tgene_id\nSample_001\tENSG00000141510\n"

    @render.download_button(filename="pathogenicity_ranks_template.tsv")
    def download_template_patho():
        yield "gene_id\tgene_symbol\tSample_001\nENSG00000141510\tTP53\t0.98\n"

    @render.ui
    def progress_bar_ui():
        val = progress_percent()
        detail = progress_detail()
        return ui.div(
            ui.div(
                ui.div(
                    f"%{val}",
                    class_="progress-bar progress-bar-striped progress-bar-animated bg-success",
                    role="progressbar",
                    style=f"width: {val}%;",
                    aria_valuenow=val,
                    aria_valuemin=0,
                    aria_valuemax=100
                ),
                class_="progress",
                style="height: 25px; font-size: 14px; font-weight: bold;"
            ),
            ui.p(detail, class_="text-muted mt-2 mb-0")
        )

    @render.code
    def pipeline_logs():
        return log_messages()

# ------------------------------------------------------------------------------
# 4. UYGULAMAYI BAŞLATMA
# ------------------------------------------------------------------------------
app = App(app_ui, server)
# iMAT & Splicing Pathogenicity Interactive Analytics Dashboard

A Python Shiny web application for interactive visualization, cohort statistical analysis, and live calculation of splicing pathogenicity and iMAT metabolic pipeline models.

---

## Overview

This repository provides an interactive analytics platform designed to visualize and process flux variability, target gene expression, and reporter-like pathway comparisons across matched sample cohorts (Control vs. Disease/AD). 

It connects directly with underlying analysis workflows and synthetic pipeline datasets from [Bioinformatikci/bioinformatics-portfolio-projects](https://github.com/Bioinformatikci/bioinformatics-portfolio-projects).

---

## Key Features

  - **Live Computed Analytics:** Dynamic evaluation of iMAT flux thresholds and corrected p-value limits.
  - **Interactive Visualizations (Plotly):**
  - **Reporter-like Pathway Comparison:** Horizontal comparative bar plots with threshold lines.
  - **Staged iMAT Reaction Activity:** Multi-stage boxplot distribution per sample group.
  - **Paired Model-Size Changes:** Direct model shrinkage/expansion analysis.
  - **Target Expression Categories:** Low, Neutral, and High target expression categorization.
  - **Subsystem-Target Network Map:** Interactive structural layout linking subsystems, genes, and samples.
  - **Data Table & Reporting:** Real-time search/filter capabilities with single-click CSV and Executive HTML report downloads.
- **Custom Dark Interface:** Sleek, high-contrast dark theme optimized for biological matrix interpretation.

---

## Project Structure

```text
├── app.py                   # Main Python Shiny Application
├── external_repo/           # Submodule / Forked core pipeline repository
│   └── bioinformatics-portfolio-projects/
└── README.md                # Project documentation

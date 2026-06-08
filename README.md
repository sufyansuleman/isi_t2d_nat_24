# ISI-T2D Association Browser

**Interactive data browser for T2D genetic variants and thier associaiton with insulin sensitivity indices:**
> Suzuki et al. *Genetic drivers of heterogeneity in type 2 diabetes pathophysiology.* **Nature** (2024).

**Live browser → https://sufyansuleman.github.io/isi_t2d_nat_24/**

---

## What this is

Genetic variants identified in the largest T2D GWAS to date (>2.5 million individuals, Suzuki et al. *Nature* 2024) were tested for association with **23 insulin sensitivity indices (ISIs)** derived from fasting blood measurements and oral glucose tolerance tests (OGTTs).

This browser lets you interactively explore:

- Which ISIs a given T2D risk allele associates with, and in which direction
- Whether effects operate through fasting insulin resistance, post-load glucose handling, or combined pathways
- How much of each genetic effect is mediated through BMI
- Side-by-side comparison of BMI-adjusted vs BMI-unadjusted effect sizes

All analyses are presented for both **BMI-adjusted** and **BMI-unadjusted** models.

---

## Browser tabs

| Tab | Description |
|-----|-------------|
| **Heatmap** | Effect-size heatmap across all variants × 23 ISIs with searchable data table |
| **BMI Comparison** | Side-by-side BMI-adjusted vs unadjusted heatmaps; optional per-variant Δβ calculator |
| **Volcano** | Volcano plots (β vs −log₁₀ p) per IS index, individually or all 23 at once |
| **Effect Profile** | Lollipop chart for a selected variant across all 23 ISIs (both analyses overlaid) |
| **BMI Mediation** | Baron–Kenny proportion-mediated bar chart (top 10 attenuated, top 10 amplified) |
| **Data Browser** | Searchable, downloadable table of all association results |
| **Gene Summary** | Aggregated heatmap and table by gene |

All views respond to the shared **search / filter panel** (gene name, rsID, index group, p-value threshold, minimum \|β\|).

---

## Methods

**Insulin sensitivity indices** span three time-point groups:

- *Fasting* — inv-FIns, inv-HOMA-IR, QUICKI, Raynaud SI, Belfiore basal, ISI basal, Bennett SI, Avignon SI₀, FIRI, inv-FIns/FGlu
- *OGTT 0–120 min* — inv-Ins 120, inv-Glu 120, ISI 120, inv-Ins/Glu120, Gutt Index, Avignon SI120, Avignon SIM
- *OGTT 0–30–120 min* — Stumvoll Modi, Stumvoll Dem, inv-IFC, BIGTT SI, Matsuda, Matsuda AUC

**BMI mediation** uses the Baron–Kenny formula:
PM = (β_total − β_direct) / |β_total| × 100, where β_total is the unadjusted effect and β_direct is the BMI-adjusted effect. Only variants with |β_total| ≥ 0.005 across ≥ 3 ISIs are included; individual PMs are clamped to ±300 % before taking the median.

---

## Citations

If you use data or visualisations from this browser, please cite both papers:

> Suzuki K, Hatzikotoulas K, Southam L, et al. *Genetic drivers of heterogeneity in type 2 diabetes pathophysiology.* **Nature** 627, 347–357 (2024). https://doi.org/10.1038/s41586-024-07019-6

> Suleman S, Ängquist L, Linneberg A, Hansen T, Grarup N. *Exploring the genetic intersection between obesity-associated genetic variants and insulin sensitivity indices.* **Scientific Reports** 15, 15761 (2025). https://doi.org/10.1038/s41598-025-98507-w

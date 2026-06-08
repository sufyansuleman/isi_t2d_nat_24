# ISI-T2D Association Browser

**Interactive supplementary data browser for:**
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

## Repository structure

```
isi_t2d_nat_24/
├── index.qmd           # Main Quarto source — all interactive tabs (OJS + Observable Plot)
├── _quarto.yml         # Site config: output to docs/, theme cosmo, body-width 1400px
├── www/
│   └── custom.css      # Page layout and control panel styles
├── R/
│   ├── load_data.R     # Loads and harmonises the TSV association files
│   └── constants.R     # IS index order and groupings
├── docs/               # Rendered static site (served by GitHub Pages)
│   └── index.html
└── data/               # Association TSVs — NOT committed (see note below)
```

> **Data files** (`data/*.txt`) contain unpublished association statistics and are excluded from this repository via `.gitignore`. The rendered `docs/index.html` embeds all necessary data at build time via Quarto's `ojs_define()`.

---

## Technical stack

- **[Quarto](https://quarto.org/)** — renders `index.qmd` to static HTML in `docs/`
- **[Observable JS](https://observablehq.com/@observablehq/inputs)** + **[Observable Plot](https://observablehq.com/plot/)** — all interactivity is client-side; no server required
- **[GitHub Pages](https://pages.github.com/)** — serves `docs/` from the `main` branch
- **R** (`dplyr`, `readr`, `here`) — data loading at render time only

To rebuild the site locally:

```r
# In RStudio or an R terminal at the project root
quarto::quarto_render("index.qmd")
```

The rendered `docs/index.html` can then be previewed locally or pushed to update the live site.

---

## Citation

If you use data or visualisations from this browser, please cite:

> Suzuki K, Hatzikotoulas K, Southam L, et al. *Genetic drivers of heterogeneity in type 2 diabetes pathophysiology.* **Nature** 627, 347–357 (2024). https://doi.org/10.1038/s41586-024-07019-6

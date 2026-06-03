# R/constants.R ─────────────────────────────────────────────────────────────
# Single source of truth for all shared constants used by load_data.R,
# plots.R, and dashboard.Rmd.

# ── Phenotype groupings ───────────────────────────────────────────────────────

FASTING_PHENOTYPES <- c(
  "insu0", "raynaud", "homa_ir", "quicki", "belfiore_0",
  "ig_ratio_0", "isi_0", "bennett", "avignon_si0", "firi"
)

OGTT_0_120_PHENOTYPES <- c(
  "insu120", "pglu120", "isi_120", "ig_ratio_120", "gutt",
  "avignon_si120", "avignon_sim", "mod_stumvoll", "stumvoll_dem", "ifc"
)

OGTT_0_30_120_PHENOTYPES <- c("bigtt_si", "matsuda", "matsuda_auc")

ALL_PHENOTYPES <- c(
  FASTING_PHENOTYPES, OGTT_0_120_PHENOTYPES, OGTT_0_30_120_PHENOTYPES
)

# ── IS Index labels and ordering ──────────────────────────────────────────────

# Named vector: phenotype_code -> human-readable label
IS_INDEX_LABELS <- c(
  insu0         = "inv-FIns",
  raynaud       = "Raynaud SI",
  homa_ir       = "inv-HOMA-IR",
  quicki        = "QUICKI",
  belfiore_0    = "Belfiore basal",
  ig_ratio_0    = "inv-FIns/FGlu",
  isi_0         = "ISI basal",
  bennett       = "Bennett SI",
  avignon_si0   = "Avignon SI0",
  insu120       = "inv-Ins 120",
  pglu120       = "inv-Glu 120",
  isi_120       = "ISI 120",
  ig_ratio_120  = "inv-Ins/Glu120",
  gutt          = "Gutt Index",
  avignon_si120 = "Avignon SI120",
  avignon_sim   = "Avignon SIM",
  mod_stumvoll  = "Stumvoll Modi",
  stumvoll_dem  = "Stumvoll Dem",
  ifc           = "inv-IFC",
  bigtt_si      = "BIGTT SI",
  matsuda       = "Matsuda",
  firi          = "FIRI",
  matsuda_auc   = "Matsuda AUC"
)

# Fixed display order for x-axis (fasting → OGTT 0,120 → OGTT 0,30,120)
IS_INDEX_ORDER <- c(
  "inv-FIns",       "inv-HOMA-IR",    "Raynaud SI",     "QUICKI",
  "Belfiore basal", "inv-FIns/FGlu",  "ISI basal",      "Bennett SI",     "Avignon SI0",
  "FIRI",
  "inv-Ins 120",    "inv-Glu 120",    "ISI 120",        "inv-Ins/Glu120", "Gutt Index",
  "Avignon SI120",  "Avignon SIM",    "Stumvoll Modi",  "Stumvoll Dem",   "inv-IFC",
  "BIGTT SI",       "Matsuda",        "Matsuda AUC"
)

# ── Colours ───────────────────────────────────────────────────────────────────

INDEX_GROUP_COLORS <- c(
  "Fasting"        = "#2c7bb6",
  "OGTT,0-120"     = "#d7191c",
  "OGTT,0-30-120"  = "#1a9641"
)

# Diverging blue-white-red scale (matches original; red = higher IS)
BETA_COLORSCALE <- list(
  list(0.00, "#053061"),
  list(0.20, "#2166ac"),
  list(0.40, "#92c5de"),
  list(0.50, "#f7f7f7"),
  list(0.60, "#f4a582"),
  list(0.80, "#d6604d"),
  list(1.00, "#67001f")
)

# Purple-white-green for BMI attenuation difference map
ATTENUATION_COLORSCALE <- list(
  list(0.00, "#762a83"),   # purple  → BMI attenuates (β shrinks with adjustment)
  list(0.40, "#c2a5cf"),
  list(0.50, "#f7f7f7"),   # white   → no change
  list(0.60, "#7fbf7b"),
  list(1.00, "#1b7837")    # green   → BMI amplifies (suppressor)
)

# Significance colour tiers
SIG_COLORS <- c(
  "p < 0.001" = "#d73027",
  "p < 0.05"  = "#fc8d59",
  "NS"        = "#cccccc"
)

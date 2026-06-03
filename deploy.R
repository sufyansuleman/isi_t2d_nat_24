# deploy.R ───────────────────────────────────────────────────────────────────
# One-click deployment of the ISI-T2D Shiny dashboard to shinyapps.io.
#
# FIRST-TIME SETUP
# ─────────────────
# 1. Create a free account at https://www.shinyapps.io
# 2. Go to Account → Tokens → Show → Copy
# 3. Paste the rsconnect::setAccountInfo() call below and run it ONCE:
#
#    rsconnect::setAccountInfo(
#      name   = "YOUR_SHINYAPPS_USERNAME",
#      token  = "YOUR_TOKEN",
#      secret = "YOUR_SECRET"
#    )
#
# DEPLOYMENT
# ───────────
# Run this whole script from the project root (or use Ctrl+Shift+S in RStudio):
#   source("deploy.R")
#
# The app URL will be:
#   https://YOUR_SHINYAPPS_USERNAME.shinyapps.io/isi-t2d-browser/
# ─────────────────────────────────────────────────────────────────────────────

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  install.packages("rsconnect")
}
library(rsconnect)

# Files to bundle — include data only if they are small enough (< 100 MB)
app_files <- c(
  "dashboard.Rmd",
  "R/constants.R",
  "R/load_data.R",
  "R/plots.R",
  "www/custom.css"
)

# Optionally include data files if they exist and aren't too large
data_files <- list.files("data", pattern = "\\.txt$|\\.tsv$|\\.csv$",
                          full.names = TRUE, recursive = FALSE)

if (length(data_files) > 0) {
  total_mb <- sum(file.size(data_files)) / 1e6
  if (total_mb < 200) {
    message("Including ", length(data_files), " data files (",
            round(total_mb, 1), " MB) in deployment bundle.")
    app_files <- c(app_files, paste0("data/", basename(data_files)))
  } else {
    message("Data files total ", round(total_mb, 1), " MB — too large to bundle.\n",
            "Consider hosting data externally and loading via URL in load_data.R.")
  }
}

rsconnect::deployApp(
  appDir       = here::here(),
  appFiles     = app_files,
  appName      = "isi-t2d-browser",
  appTitle     = "ISI-T2D Genetic Association Browser",
  launch.browser = TRUE,
  forceUpdate  = TRUE
)

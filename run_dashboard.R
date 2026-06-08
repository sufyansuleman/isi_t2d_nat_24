setwd("C:/Users/rnh585/Documents/Manuscripts_Supplementray_data/isi_t2d_nat_24")

# Locate pandoc — RStudio sets RSTUDIO_PANDOC; outside RStudio we search common locations
if (!rmarkdown::pandoc_available()) {
  candidate_dirs <- c(
    Sys.getenv("RSTUDIO_PANDOC"),
    file.path(Sys.getenv("LOCALAPPDATA"), "Programs/RStudio/resources/app/bin/quarto/bin/tools"),
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools",
    file.path(Sys.getenv("PROGRAMFILES"), "RStudio/resources/app/bin/quarto/bin/tools")
  )
  for (d in candidate_dirs[nchar(candidate_dirs) > 0]) {
    if (file.exists(file.path(d, "pandoc.exe"))) {
      rmarkdown::find_pandoc(dir = d, cache = FALSE)
      break
    }
  }
}

if (!rmarkdown::pandoc_available()) {
  stop(
    "pandoc not found.\n",
    "  Option 1: Open dashboard.Rmd in RStudio and click 'Run Document'.\n",
    "  Option 2: Install pandoc from https://pandoc.org and add it to PATH."
  )
}

options(shiny.port = 3939, shiny.host = "127.0.0.1")
rmarkdown::run(
  "dashboard.Rmd",
  shiny_args = list(port = 3939, host = "127.0.0.1", launch.browser = TRUE)
)


library(here)

# create folders for the R project
dirs <- c(
  "data/processed",
  "scripts",
  "results/tables",
  "results/figures",
  "config",
  "docs"
)

for (d in dirs) {
  dir.create(here(d), recursive = TRUE, showWarnings = FALSE)
}
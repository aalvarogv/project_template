#!/usr/bin/env bash
set -euo pipefail

# --- R (librería local sin sudo) ---
mkdir -p ./R_LIBS
Rscript -e 'repos <- "https://cloud.r-project.org";
            dir.create("R_LIBS", showWarnings = FALSE);
            pkgs <- c("optparse","igraph","tidyverse","jsonlite","gprofiler2");
            inst <- pkgs[!(pkgs %in% rownames(installed.packages(lib.loc = "./R_LIBS")))];
            if (length(inst)) install.packages(inst, lib="./R_LIBS", repos=repos)'

# --- Python ---
python -m venv .venv || true
. .venv/bin/activate || . .venv/Scripts/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

echo "[ok] R y Python listos (./R_LIBS y .venv)"

# Crypto_POT_model

This repository implements an end-to-end empirical pipeline for **tail risk modeling and risk-spillover analysis** on **cryptocurrency daily price series**.

The workflow follows five stages:

- **Step 0**: Prepare a panel of daily prices.
- **Step 1**: Estimate dynamic POT-GPD tail parameters by MLE (tail index and scale), and compare POT scale dynamics with GARCH(1,1)-based scale dynamics.
- **Step 2**: Feed tail-index series into the Diebold–Yilmaz (DY) connectedness framework, run FEVD, obtain spillover matrices, and build risk-spillover networks.
- **Step 3**: Run robustness checks over DY rolling-window size and POT hurdle (threshold quantile).
- **Step 4**: Batch-convert outputs into LaTeX-ready figures/tables.

---

## 1) Repository structure

```text
Crypto_POT_model/
├── data_processing/        # Step 0: build panel_dt from raw CSV files
├── POT/                    # Step 1: dynamic POT-GPD estimation + plots
├── DY/                     # Step 2: DY connectedness / FEVD / spillover outputs
├── Robustness/             # Step 3: threshold and window-size robustness
├── Formatting/             # Step 4: batch LaTeX export helpers
├── Params/                 # save parameters
├── Tex/                    # TeX artifacts/templates
└── tools/                  # utility scripts
```

---

## 2) Requirements

The project is written in **R**.

### 2.1 R packages used across scripts

- Data & utilities: `data.table`, `zoo`, `xts`, `lubridate`, `stringr`, `dplyr`
- Statistics/econometrics: `rugarch`, `rmgarch`, `urca`, `tseries`, `TSA`, `moments`, `fBasics`
- Optimization/parallel: `GenSA`, `numDeriv`, `foreach`, `doParallel`
- Plotting/export: `ggplot2`, `scales`, `svglite`, `igraph`, `openxlsx`, `writexl`, `readxl`
- Other helper packages: `progress`, `rcompanion`

Install missing packages:

```r
packages <- c(
  "data.table", "zoo", "xts", "lubridate", "stringr", "dplyr",
  "rugarch", "rmgarch", "urca", "tseries", "TSA", "moments", "fBasics",
  "GenSA", "numDeriv", "foreach", "doParallel",
  "ggplot2", "scales", "svglite", "igraph", "openxlsx", "writexl", "readxl",
  "progress", "rcompanion"
)

missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing)
```

> Note: `DY/main.R` calls `ConnectednessApproach(...)`. Make sure the corresponding connectedness package/function is installed and available in your R environment.

---

## 3) Step 0 — Data preparation

**Script**: `data_processing/data_processing.R`

### Input

- Raw CSV files under `data/ori_data/`
- Each file should include at least: `End`, `Close`

### Processing

- Read each coin CSV file.
- Parse `End` as date and `Close` as numeric.
- Rename each `Close` column to coin name.
- Merge all series by date into one wide panel.

### Output

- `data/panel_dt.RData`

---

## 4) Step 1 — Dynamic POT-GPD estimation (MLE)

**Main script**: `POT/main.R`

### What it does

- Computes log returns from prices (scaled by 100).
- Builds downside/upside exceedance series relative to quantile hurdle.
- Estimates dynamic POT-GPD parameters by MLE:
  - `xi` (tail-index parameter in the script)
  - `sigma` (scale parameter)
- Saves diagnostics and plots per coin.

### Key configurable options (top of script)

- `q = 0.95` (baseline hurdle; robustness examples: `0.90`, `0.85`)
- `target_folder` output root
- `neg.only` / `pos.only` to run one tail only
- parameter-bound presets loaded from `Params/est_param_*.R`

### Typical outputs

- `result/POT/percent95/ret_desc_stats/*.RData`
- `result/POT/percent95/down_tail/optim_result/*.RData`
- `result/POT/percent95/up_tail/optim_result/*.RData`
- Plot PDFs (`P_xi_exc`, `P_std_sigma`, `return_plot`, `price_plot`)

### POT scale vs GARCH(1,1) scale comparison

The project already depends on `rugarch`, so you can compare POT `sigma` against GARCH(1,1) conditional volatility (or scale proxy) on aligned samples. Common checks:

- summary-stat comparison (mean/median/quantiles)
- time-series correlation/co-movement
- behavior during stress windows (extreme-return days)

---

## 5) Step 2 — DY connectedness, FEVD, and spillovers

**Main script**: `DY/main.R`

### What it does

1. Loads Step 1 per-coin estimation objects (up/down tails).
2. Extracts tail-index series and transforms to `1/xi` for connectedness input.
3. Runs DY connectedness/FEVD computations.
4. Produces spillover objects/matrices (`dca_up`, `dca_down`) for downstream analysis.

### Key configurable options

- `window.size = 200` (robustness: `150`, `250`)
- `model = "VAR"`
- `nlag = 1`, `nfore = 20`
- `core.mode = FALSE` (set `TRUE` for a core-coin subset)

### Output

- `result/DY/percent95win200/dca.RData` (or `dca_core.RData`)

Useful plotting/analysis scripts:

- `DY/plot.R`
- `DY/Risk_spillovers_plot.R`

---

## 6) Step 3 — Robustness checks

**Main script**: `Robustness/main.R`

### Dimensions tested

1. **POT hurdle robustness**: `0.95` (baseline) vs `0.90` vs `0.85`
2. **DY window robustness**: `200` (baseline) vs `150` vs `250`

### Typical outputs

- `result/Robustness/TCI_plot/*.pdf` (total connectedness index comparisons)
- `result/Robustness/rank45plot/*.pdf` (net-spillover rank comparison plots)

---

## 7) Step 4 — Batch LaTeX export

- Figure export helper: `Formatting/figure.R`
- Table export helper: `Formatting/table.R`

### What gets generated

- LaTeX-ready image blocks (e.g., minipage grids)
- LaTeX tables from saved R objects (descriptive stats, parameter estimates, spillovers)

### Example artifacts

- `table_txt/tail_index_images.txt`
- `table_txt/up_NET_plot.txt`
- `table_txt/down_NET_plot.txt`
- `table_txt/desc_table.txt`
- `table_txt/up_optim_result_table.txt`
- `table_txt/down_optim_result_table.txt`
- `table_txt/spillover_down_table.txt`
- `table_txt/spillover_up_table.txt`

---

## 8) Recommended execution order

> Important: several scripts currently use hard-coded paths such as `setwd("D:/MyFiles/EVT/Crypto_POT")`. Update these to your local path (or refactor to relative paths) before running.

1. `data_processing/data_processing.R` (skip if `panel_dt.RData` already exists)
2. `POT/main.R` (baseline `q = 0.95`)
3. `DY/main.R` (baseline `window.size = 200`)
4. `Robustness/main.R` (after robustness runs are available)
5. `Formatting/table.R` and `Formatting/figure.R`

---

## 9) Reading the outputs

- **Tail index / `1/xi`**: reflects tail heaviness and extreme-risk exposure.
- **NET spillover**: positive = net transmitter; negative = net receiver.
- **TCI**: aggregate market-wide spillover intensity.
- **Robustness consistency** across thresholds/windows strengthens inference credibility.

---

## 10) Reproducibility and extension notes

- Keep `set.seed(...)` fixed and log software/package versions.
- Add run logs/error traces for long multi-coin loops.
- Potential extensions:
  - regime-split estimation (bull/bear or high/low volatility)
  - exogenous controls (macro factors, on-chain metrics)
  - model benchmarking against alternative EVT approaches


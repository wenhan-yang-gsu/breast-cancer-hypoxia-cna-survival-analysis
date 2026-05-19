## Analysis code for TCGA and METABRIC hypoxia, survival, and CNA analyses.
## This script reads the cBioPortal-format files, prepares the analysis datasets,
## runs the survival and CNA analyses, and exports the main tables and figures.
##
## Reproducibility notes:
##   This script expects a project directory containing:
##     brca_tcga_pan_can_atlas_2018/
##     brca_metabric/
##
##   You may set `project_dir` manually below.
##   If left as NULL, the script will try to find the project directory automatically.
##
##   Running this script will create:
##     outputs_final/
##
## Software environment:
##   The analysis was run in R 4.5.2. See sessionInfo.txt for details.
## =====================================================================

options(stringsAsFactors = FALSE)

## ---------------------------------------------------------------------
## 0) User configuration
## ---------------------------------------------------------------------

## Option A: set manually if you want
project_dir <- NULL
## Example:
## project_dir <- "/Users/yourname/Desktop/Dataset"

## Option B: auto-detect if project_dir is NULL
find_project_dir <- function() {
  candidates <- c(
    getwd(),
    path.expand("~/Desktop/Dataset"),
    dirname(getwd())
  )
  
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
  
  for (d in candidates) {
    if (dir.exists(file.path(d, "brca_tcga_pan_can_atlas_2018")) &&
        dir.exists(file.path(d, "brca_metabric"))) {
      return(d)
    }
  }
  
  NA_character_
}

if (is.null(project_dir)) {
  project_dir <- find_project_dir()
}

if (is.na(project_dir) || !dir.exists(project_dir)) {
  stop(
    paste0(
      "Cannot find the project directory automatically.\n",
      "Please set `project_dir` manually to the folder containing both:\n",
      "  brca_tcga_pan_can_atlas_2018/\n",
      "  brca_metabric/\n"
    )
  )
}

setwd(project_dir)
cat("Working directory:", getwd(), "\n")

tcga_dir     <- "brca_tcga_pan_can_atlas_2018"
metabric_dir <- "brca_metabric"

if (!dir.exists(tcga_dir))     stop("TCGA folder not found: ", tcga_dir)
if (!dir.exists(metabric_dir)) stop("METABRIC folder not found: ", metabric_dir)

out_dir   <- file.path(getwd(), "outputs_final")
plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
data_dir  <- file.path(out_dir, "processed_data")
diag_dir  <- file.path(out_dir, "diagnostics")

dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir,  recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------
## 1) Packages
## ---------------------------------------------------------------------

req_pkgs <- c(
  "data.table", "dplyr", "stringr", "tidyr", "purrr", "readr",
  "survival", "survminer", "ggplot2", "broom", "openxlsx", "tibble"
)

to_install <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

invisible(lapply(req_pkgs, library, character.only = TRUE))

theme_set(theme_bw(base_size = 12))

model_registry <- list()
ph_registry    <- list()
skip_registry  <- list()

## ---------------------------------------------------------------------
## 2) Helper functions
## ---------------------------------------------------------------------

read_cbio_table <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  x <- readLines(path, warn = FALSE)
  idx <- which(!grepl("^#", x) & grepl("\t", x))[1]
  data.table::fread(path, skip = idx - 1, sep = "\t", data.table = FALSE)
}

pick_first_file <- function(dir, candidates) {
  for (fn in candidates) {
    fp <- file.path(dir, fn)
    if (file.exists(fp)) return(fp)
  }
  NA_character_
}

pick_col <- function(df, patterns) {
  nms <- names(df)
  for (p in patterns) {
    hit <- nms[grepl(p, nms, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

pick_all_cols <- function(df, patterns) {
  nms <- names(df)
  unique(unlist(lapply(patterns, function(p) nms[grepl(p, nms, ignore.case = TRUE)])))
}

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

status_to_event <- function(x) {
  x0 <- trimws(toupper(as.character(x)))
  out <- rep(NA_integer_, length(x0))
  
  suppressWarnings(num <- as.numeric(x0))
  out[!is.na(num)] <- ifelse(num[!is.na(num)] == 1, 1L, 0L)
  
  out[is.na(out) & grepl("DECEASED|DEAD|DIED|EVENT", x0)] <- 1L
  out[is.na(out) & grepl("LIVING|ALIVE|CENSORED|NO EVENT", x0)] <- 0L
  out[is.na(out) & grepl("^1[:_]", x0)] <- 1L
  out[is.na(out) & grepl("^0[:_]", x0)] <- 0L
  
  out
}

pick_survival_cols <- function(df) {
  t_os <- pick_col(df, c("^OS_MONTHS$", "^OS_TIME$", "^OS_MONTH$", "OS_MONTH"))
  s_os <- pick_col(df, c("^OS_STATUS$", "^OS_EVENT$", "^OS_IND$", "OS_STATUS"))
  if (!is.na(t_os) && !is.na(s_os)) return(list(time = t_os, status = s_os, endpoint = "OS"))
  
  t_dss <- pick_col(df, c("^DSS_MONTHS$", "^DSS_TIME$", "^DSS_MONTH$"))
  s_dss <- pick_col(df, c("^DSS_STATUS$", "DSS_STATUS"))
  if (!is.na(t_dss) && !is.na(s_dss)) return(list(time = t_dss, status = s_dss, endpoint = "DSS"))
  
  t_dfs <- pick_col(df, c("^DFS_MONTHS$", "^DFS_TIME$", "^DFS_MONTH$"))
  s_dfs <- pick_col(df, c("^DFS_STATUS$", "DFS_STATUS"))
  if (!is.na(t_dfs) && !is.na(s_dfs)) return(list(time = t_dfs, status = s_dfs, endpoint = "DFS"))
  
  t_pfs <- pick_col(df, c("^PFS_MONTHS$", "^PFS_TIME$", "^PFS_MONTH$"))
  s_pfs <- pick_col(df, c("^PFS_STATUS$", "PFS_STATUS"))
  if (!is.na(t_pfs) && !is.na(s_pfs)) return(list(time = t_pfs, status = s_pfs, endpoint = "PFS"))
  
  stop("No survival endpoint columns found.")
}

pick_buffa_col <- function(df) {
  pick_col(df, c("BUFFA", "HYPOXIA.*BUFFA", "BUFFA.*HYPOXIA"))
}

pick_subtype_col_tcga <- function(df) {
  pick_col(df, c("^SUBTYPE$", "PAM50", "INTRINSIC_SUBTYPE", "CLAUDIN", "BREAST_SUBTYPE"))
}

pick_subtype_col_metabric <- function(df) {
  pick_col(df, c("^CLAUDIN_SUBTYPE$", "PAM50", "INTRINSIC", "^SUBTYPE$", "CLAUDIN", "HISTOLOGICAL_SUBTYPE"))
}

pick_age_col <- function(df) {
  pick_col(df, c("^AGE$", "AGE_AT", "PATIENT_AGE", "AGE_AT_DIAGNOSIS", "AGE_AT_INITIAL_PATHOLOGIC_DIAGNOSIS"))
}

pick_stage_col <- function(df) {
  pick_col(df, c("AJCC.*STAGE", "PATHOLOGIC.*STAGE", "TUMOR_STAGE", "^STAGE$", "CLINICAL_STAGE"))
}

pick_grade_col <- function(df) {
  pick_col(df, c("GRADE", "HISTOLOGIC.*GRADE", "TUMOR_GRADE", "^GRADE$"))
}

simplify_stage <- function(stage_chr) {
  x <- trimws(toupper(as.character(stage_chr)))
  suppressWarnings(num <- as.numeric(x))
  
  out <- dplyr::case_when(
    !is.na(num) & num == 4 ~ "IV",
    !is.na(num) & num == 3 ~ "III",
    !is.na(num) & num == 2 ~ "II",
    !is.na(num) & num == 1 ~ "I",
    !is.na(num) & num == 0 ~ NA_character_,
    grepl("\\bIV\\b", x) ~ "IV",
    grepl("\\bIII\\b", x) ~ "III",
    grepl("\\bII\\b", x) ~ "II",
    grepl("(^|[^V])I([^I]|$)", x) ~ "I",
    TRUE ~ NA_character_
  )
  
  factor(out, levels = c("I", "II", "III", "IV"))
}

stage_to_2cat <- function(stage_simple) {
  x <- as.character(stage_simple)
  y <- dplyr::case_when(
    x %in% c("I", "II") ~ "I_II",
    x %in% c("III", "IV") ~ "III_IV",
    TRUE ~ NA_character_
  )
  factor(y, levels = c("I_II", "III_IV"))
}

simplify_grade <- function(grade_chr) {
  x <- trimws(toupper(as.character(grade_chr)))
  suppressWarnings(num <- as.numeric(x))
  
  out <- dplyr::case_when(
    !is.na(num) & num == 1 ~ "G1",
    !is.na(num) & num == 2 ~ "G2",
    !is.na(num) & num == 3 ~ "G3",
    grepl("GRADE ?1|\\bG1\\b|LOW", x) ~ "G1",
    grepl("GRADE ?2|\\bG2\\b|INTERMEDIATE|MODERATE", x) ~ "G2",
    grepl("GRADE ?3|\\bG3\\b|HIGH|POOR", x) ~ "G3",
    TRUE ~ NA_character_
  )
  
  factor(out, levels = c("G1", "G2", "G3"))
}

grade_to_2cat <- function(grade_simple) {
  x <- as.character(grade_simple)
  y <- dplyr::case_when(
    x %in% c("G1", "G2") ~ "G1_G2",
    x %in% c("G3") ~ "G3",
    TRUE ~ NA_character_
  )
  factor(y, levels = c("G1_G2", "G3"))
}

text_to_yes_no <- function(x) {
  x0 <- trimws(toupper(as.character(x)))
  out <- rep(NA_character_, length(x0))
  
  suppressWarnings(num <- as.numeric(x0))
  out[!is.na(num) & num == 1] <- "Yes"
  out[!is.na(num) & num == 0] <- "No"
  
  out[grepl("^YES$|RECEIVED|ADMINISTERED|TREATED|COMPLETED|PERFORMED|DONE", x0)] <- "Yes"
  out[grepl("^NO$|NONE|NOT ADMINISTERED|NOT RECEIVED|UNTREATED|NAIVE|NEVER", x0)] <- "No"
  out[grepl("UNKNOWN|NOT AVAILABLE|N/A|NA|MISSING|UNSPECIFIED", x0)] <- NA_character_
  
  factor(out, levels = c("No", "Yes"))
}

pick_treatment_cols <- function(df) {
  pick_all_cols(df, c("TREAT", "THERAP", "CHEMO", "RADIAT", "HORMON", "ENDOCRINE", "TARGET"))
}

derive_any_treatment <- function(df) {
  cols <- pick_treatment_cols(df)
  if (length(cols) == 0) {
    return(factor(rep(NA_character_, nrow(df)), levels = c("No", "Yes")))
  }
  
  yn_list <- lapply(cols, function(cc) as.character(text_to_yes_no(df[[cc]])))
  any_obs <- Reduce(`|`, lapply(yn_list, function(v) !is.na(v)))
  any_yes <- Reduce(`|`, lapply(yn_list, function(v) !is.na(v) & v == "Yes"))
  
  out <- rep(NA_character_, nrow(df))
  out[any_obs & !any_yes] <- "No"
  out[any_yes] <- "Yes"
  
  factor(out, levels = c("No", "Yes"))
}

count_complete_model <- function(df, vars) {
  vars <- unique(c("time", "event", vars))
  vars <- vars[vars %in% names(df)]
  d <- df %>% dplyr::select(all_of(vars)) %>% filter(complete.cases(.), time > 0)
  data.frame(
    n_complete = nrow(d),
    events_complete = sum(d$event == 1, na.rm = TRUE)
  )
}

covariate_availability <- function(df, vars, file_stub = NULL) {
  out <- lapply(vars, function(v) {
    if (!v %in% names(df)) {
      return(data.frame(
        variable = v,
        exists = FALSE,
        n_nonmissing = 0,
        prop_complete = 0,
        n_unique = 0,
        usable = FALSE
      ))
    }
    
    x <- df[[v]]
    data.frame(
      variable = v,
      exists = TRUE,
      n_nonmissing = sum(!is.na(x)),
      prop_complete = round(mean(!is.na(x)), 3),
      n_unique = dplyr::n_distinct(x[!is.na(x)]),
      usable = mean(!is.na(x)) > 0 && dplyr::n_distinct(x[!is.na(x)]) >= 2
    )
  }) %>% bind_rows()
  
  if (!is.null(file_stub)) {
    readr::write_csv(out, file.path(diag_dir, paste0(file_stub, "_covariate_availability.csv")))
  }
  
  out
}

available_covars <- function(df, vars, min_complete_prop = 0.40) {
  tab <- covariate_availability(df, vars)
  tab %>%
    filter(exists, usable, prop_complete >= min_complete_prop) %>%
    pull(variable)
}

save_skip_note <- function(model_name, reason, out_dir) {
  note <- data.frame(model = model_name, reason = reason)
  readr::write_csv(note, file.path(out_dir, paste0(model_name, "_SKIPPED.csv")))
  skip_registry[[length(skip_registry) + 1L]] <<- note
}

save_model_bundle <- function(fit, data_used, model_name, out_dir, formula_text = NULL) {
  model_info <- data.frame(
    model = model_name,
    n = nrow(data_used),
    events = sum(data_used$event == 1, na.rm = TRUE),
    formula = if (is.null(formula_text)) deparse(formula(fit)) else formula_text
  )
  
  readr::write_csv(model_info, file.path(out_dir, paste0(model_name, "_model_info.csv")))
  
  tidy_tab <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  names(tidy_tab)[names(tidy_tab) == "estimate"]  <- "HR"
  names(tidy_tab)[names(tidy_tab) == "conf.low"]  <- "CI_low"
  names(tidy_tab)[names(tidy_tab) == "conf.high"] <- "CI_high"
  readr::write_csv(tidy_tab, file.path(out_dir, paste0(model_name, "_tidy.csv")))
  
  capture.output(summary(fit), file = file.path(out_dir, paste0(model_name, "_summary.txt")))
  
  ph <- survival::cox.zph(fit)
  ph_df <- as.data.frame(ph$table)
  ph_df$term  <- rownames(ph_df)
  rownames(ph_df) <- NULL
  ph_df$model <- model_name
  readr::write_csv(ph_df, file.path(out_dir, paste0(model_name, "_PH_test.csv")))
  
  pdf(file.path(out_dir, paste0(model_name, "_PH_plots.pdf")), width = 8, height = 6)
  plot(ph)
  dev.off()
  
  model_registry[[length(model_registry) + 1L]] <<- model_info
  ph_registry[[length(ph_registry) + 1L]] <<- ph_df
  
  invisible(list(info = model_info, tidy = tidy_tab, ph = ph_df))
}

run_cox_model <- function(df, rhs_formula, required_vars, model_name, out_dir,
                          relevel_list = list(), min_n = 30, min_events = 10) {
  keep_vars <- unique(c("time", "event", required_vars))
  keep_vars <- keep_vars[keep_vars %in% names(df)]
  
  d <- df %>%
    dplyr::select(all_of(keep_vars)) %>%
    filter(complete.cases(.), time > 0)
  
  if (nrow(d) < min_n) {
    save_skip_note(model_name, paste0("Too few complete cases: n = ", nrow(d), " < ", min_n), out_dir)
    return(NULL)
  }
  
  if (sum(d$event == 1, na.rm = TRUE) < min_events) {
    save_skip_note(model_name, paste0("Too few events: ", sum(d$event == 1, na.rm = TRUE), " < ", min_events), out_dir)
    return(NULL)
  }
  
  for (nm in names(relevel_list)) {
    if (nm %in% names(d)) {
      d[[nm]] <- factor(d[[nm]])
      ref <- relevel_list[[nm]]
      if (ref %in% levels(d[[nm]])) {
        d[[nm]] <- relevel(d[[nm]], ref = ref)
      }
    }
  }
  
  formula_text <- paste0("Surv(time, event) ~ ", rhs_formula)
  fit <- coxph(as.formula(formula_text), data = d, x = TRUE, y = TRUE)
  save_model_bundle(fit, d, model_name, out_dir, formula_text = formula_text)
  
  invisible(list(fit = fit, data = d, formula = formula_text))
}

save_model_plan <- function(df, model_name, terms, out_dir) {
  cc <- count_complete_model(df, terms)
  out <- data.frame(
    model = model_name,
    rhs_terms = paste(terms, collapse = " + "),
    n_complete = cc$n_complete,
    events_complete = cc$events_complete
  )
  readr::write_csv(out, file.path(out_dir, paste0(model_name, "_planned_complete_cases.csv")))
  out
}

save_km_plot <- function(df, formula_text, file_stub, legend_title = "Group") {
  frm <- as.formula(formula_text, env = parent.frame())
  fit <- survival::survfit(frm, data = df)
  fit$call$formula <- frm
  
  p <- survminer::ggsurvplot(
    fit = fit,
    data = df,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    legend.title = legend_title
  )
  
  ggsave(file.path(plot_dir, paste0(file_stub, ".png")), p$plot, width = 7, height = 6, dpi = 300)
  ggsave(file.path(plot_dir, paste0(file_stub, "_risk_table.png")), p$table, width = 7, height = 3, dpi = 300)
  invisible(p)
}

calc_cna_burden_continuous <- function(cna_df) {
  if (is.null(cna_df)) return(NULL)
  
  drop_cols <- intersect(c("Hugo_Symbol", "Entrez_Gene_Id"), names(cna_df))
  sample_cols <- setdiff(names(cna_df), drop_cols)
  
  mat <- cna_df[, sample_cols, drop = FALSE]
  mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(as.character(x))))
  mat <- as.matrix(mat)
  
  data.frame(
    SAMPLE_ID = sample_cols,
    CNA_mean_abs = colMeans(abs(mat), na.rm = TRUE),
    CNA_prop_02 = colMeans(abs(mat) >= 0.2, na.rm = TRUE)
  )
}

calc_cna_burden_discrete <- function(cna_df) {
  if (is.null(cna_df)) return(NULL)
  
  drop_cols <- intersect(c("Hugo_Symbol", "Entrez_Gene_Id"), names(cna_df))
  sample_cols <- setdiff(names(cna_df), drop_cols)
  
  mat <- cna_df[, sample_cols, drop = FALSE]
  mat[] <- lapply(mat, function(x) suppressWarnings(as.numeric(as.character(x))))
  mat <- as.matrix(mat)
  
  data.frame(
    SAMPLE_ID = sample_cols,
    CNA_mean_abs = colMeans(abs(mat), na.rm = TRUE),
    CNA_prop_nonzero = colMeans(mat != 0, na.rm = TRUE)
  )
}

get_nonsyn_mut <- function(mut_df, sample_ids) {
  if (is.null(mut_df)) return(NULL)
  
  sample_col <- intersect(c("Tumor_Sample_Barcode", "SAMPLE_ID", "Sample_ID", "sample_id"), names(mut_df))[1]
  gene_col   <- intersect(c("Hugo_Symbol", "HUGO_SYMBOL", "Gene", "gene"), names(mut_df))[1]
  vc_col     <- intersect(c("Variant_Classification", "VARIANT_CLASSIFICATION"), names(mut_df))[1]
  
  if (is.na(sample_col) || is.na(gene_col)) stop("Cannot identify mutation columns.")
  
  non_syn <- c(
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "Splice_Site", "In_Frame_Del", "In_Frame_Ins", "Nonstop_Mutation",
    "Translation_Start_Site"
  )
  
  x <- mut_df %>%
    filter(.data[[sample_col]] %in% sample_ids) %>%
    transmute(
      SAMPLE_ID = as.character(.data[[sample_col]]),
      gene = as.character(.data[[gene_col]]),
      VC = if (!is.na(vc_col)) as.character(.data[[vc_col]]) else NA_character_
    ) %>%
    filter(!is.na(gene), gene != "")
  
  if (!is.na(vc_col)) {
    x <- x %>% filter(VC %in% non_syn)
  }
  
  x %>% distinct(SAMPLE_ID, gene)
}

fisher_mutOR_high <- function(tab_high_low_by01) {
  tab <- tab_high_low_by01[, c("0", "1"), drop = FALSE]
  
  a <- tab["High", "1"]; b <- tab["High", "0"]
  c <- tab["Low",  "1"]; d <- tab["Low",  "0"]
  
  a2 <- a + 0.5; b2 <- b + 0.5; c2 <- c + 0.5; d2 <- d + 0.5
  
  mutOR <- (a2 / b2) / (c2 / d2)
  logOR <- log(mutOR)
  se <- sqrt(1 / a2 + 1 / b2 + 1 / c2 + 1 / d2)
  ci <- exp(logOR + c(-1, 1) * 1.96 * se)
  
  ft <- fisher.test(tab)
  list(p = ft$p.value, mutOR = mutOR, L95 = ci[1], U95 = ci[2])
}

pick_gene_col <- function(df) {
  pick_col(df, c("^HUGO_SYMBOL$", "^Hugo_Symbol$", "^GENE$", "^Gene$"))
}

## ---------------------------------------------------------------------
## 3) Load TCGA
## ---------------------------------------------------------------------

cat("\n=============================\nLoading TCGA\n=============================\n")

tcga_patient_path <- pick_first_file(tcga_dir, c("data_clinical_patient.txt"))
tcga_sample_path  <- pick_first_file(tcga_dir, c("data_clinical_sample.txt"))
tcga_hypoxia_path <- pick_first_file(tcga_dir, c("data_clinical_supp_hypoxia.txt", "data_clinical_supp_hypoxia_score.txt"))
tcga_mut_path     <- pick_first_file(tcga_dir, c("data_mutations_extended.txt", "data_mutations.txt"))
tcga_cna_path     <- pick_first_file(tcga_dir, c("data_log2_cna.txt"))

tcga_patient <- read_cbio_table(tcga_patient_path)
tcga_sample  <- read_cbio_table(tcga_sample_path)
tcga_hypoxia <- read_cbio_table(tcga_hypoxia_path)
tcga_mut     <- if (!is.na(tcga_mut_path)) read_cbio_table(tcga_mut_path) else NULL
tcga_cna     <- if (!is.na(tcga_cna_path)) read_cbio_table(tcga_cna_path) else NULL

tcga_master <- tcga_sample %>%
  left_join(tcga_patient, by = "PATIENT_ID") %>%
  left_join(tcga_hypoxia, by = "PATIENT_ID")

tcga_surv_cols   <- pick_survival_cols(tcga_master)
tcga_buffa_col   <- pick_buffa_col(tcga_master)
tcga_subtype_col <- pick_subtype_col_tcga(tcga_master)
tcga_age_col     <- pick_age_col(tcga_master)
tcga_stage_col   <- pick_stage_col(tcga_master)
tcga_grade_col   <- pick_grade_col(tcga_master)

if (is.na(tcga_buffa_col)) stop("TCGA BUFFA hypoxia score column not found.")

tcga_df <- tcga_master %>%
  transmute(
    PATIENT_ID = as.character(PATIENT_ID),
    SAMPLE_ID  = as.character(SAMPLE_ID),
    time       = to_num(.data[[tcga_surv_cols$time]]),
    event      = status_to_event(.data[[tcga_surv_cols$status]]),
    BUFFA      = to_num(.data[[tcga_buffa_col]]),
    SUBTYPE    = if (!is.na(tcga_subtype_col)) as.character(.data[[tcga_subtype_col]]) else NA_character_,
    AGE        = if (!is.na(tcga_age_col)) to_num(.data[[tcga_age_col]]) else NA_real_,
    STAGE_RAW  = if (!is.na(tcga_stage_col)) as.character(.data[[tcga_stage_col]]) else NA_character_,
    GRADE_RAW  = if (!is.na(tcga_grade_col)) as.character(.data[[tcga_grade_col]]) else NA_character_
  ) %>%
  mutate(
    STAGE_SIMPLE = simplify_stage(STAGE_RAW),
    STAGE_2      = stage_to_2cat(STAGE_SIMPLE),
    GRADE_SIMPLE = simplify_grade(GRADE_RAW),
    GRADE_2      = grade_to_2cat(GRADE_SIMPLE)
  )

tcga_df$ANY_TREATMENT <- derive_any_treatment(tcga_master)

tcga_df <- tcga_df %>%
  filter(!is.na(time), time > 0, !is.na(event), !is.na(BUFFA)) %>%
  mutate(
    SUBTYPE = na_if(SUBTYPE, ""),
    BUFFA_z = as.numeric(scale(BUFFA)),
    hypoxia_group_glb = ifelse(BUFFA >= median(BUFFA, na.rm = TRUE), "High", "Low"),
    hypoxia_group_glb = factor(hypoxia_group_glb, levels = c("High", "Low"))
  ) %>%
  group_by(SUBTYPE) %>%
  mutate(
    hypoxia_group_st = ifelse(BUFFA >= median(BUFFA, na.rm = TRUE), "High", "Low")
  ) %>%
  ungroup() %>%
  mutate(hypoxia_group_st = factor(hypoxia_group_st, levels = c("High", "Low")))

readr::write_csv(tcga_df, file.path(data_dir, "TCGA_processed_survival_dataset.csv"))

cat("TCGA rows:", nrow(tcga_df), "\n")
cat("TCGA events:", sum(tcga_df$event == 1, na.rm = TRUE), "\n")
cat("TCGA endpoint:", tcga_surv_cols$endpoint, "\n")
cat("TCGA subtype counts:\n")
print(table(tcga_df$SUBTYPE, useNA = "ifany"))

covariate_availability(
  tcga_df,
  vars = c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"),
  file_stub = "TCGA_overall"
)

## ---------------------------------------------------------------------
## 4) TCGA global analyses
## ---------------------------------------------------------------------

cat("\n=============================\nTCGA global analyses\n=============================\n")

tcga_subtyped <- tcga_df %>%
  filter(!is.na(SUBTYPE)) %>%
  mutate(SUBTYPE = factor(SUBTYPE))

save_km_plot(
  df = tcga_df,
  formula_text = "Surv(time, event) ~ hypoxia_group_glb",
  file_stub = "TCGA_KM_global_median",
  legend_title = "TCGA hypoxia"
)

tab_glb <- table(tcga_subtyped$hypoxia_group_glb, tcga_subtyped$SUBTYPE, useNA = "no")
chi_glb <- chisq.test(tab_glb)

imbalance_df <- as.data.frame.matrix(tab_glb)
imbalance_df$Group <- rownames(imbalance_df)
rownames(imbalance_df) <- NULL
readr::write_csv(imbalance_df, file.path(table_dir, "TCGA_global_hypoxia_by_subtype_counts.csv"))
capture.output(chi_glb, file = file.path(table_dir, "TCGA_global_hypoxia_by_subtype_chisq.txt"))

prop_glb <- prop.table(tab_glb, margin = 2)
prop_glb_df <- as.data.frame(as.table(prop_glb))
names(prop_glb_df) <- c("HypoxiaGroup", "Subtype", "ProportionWithinSubtype")

p_bar <- ggplot(prop_glb_df, aes(x = Subtype, y = ProportionWithinSubtype, fill = HypoxiaGroup)) +
  geom_col(position = "stack") +
  labs(
    title = "TCGA: distribution of global-median hypoxia groups across subtypes",
    y = "Proportion within subtype", x = NULL
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(plot_dir, "TCGA_global_hypoxia_by_subtype_barplot.png"), p_bar, width = 8, height = 5, dpi = 300)

save_model_plan(tcga_subtyped, "TCGA_pooled_subtype_adjusted", c("BUFFA_z", "SUBTYPE"), diag_dir)
m_tcga_subtype <- run_cox_model(
  df = tcga_subtyped,
  rhs_formula = "BUFFA_z + SUBTYPE",
  required_vars = c("BUFFA_z", "SUBTYPE"),
  model_name = "TCGA_pooled_subtype_adjusted",
  out_dir = table_dir
)

tcga_pool_avail <- covariate_availability(
  tcga_subtyped,
  vars = c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"),
  file_stub = "TCGA_pooled"
)

tcga_pool_covars <- intersect(
  c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"),
  available_covars(tcga_subtyped, c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"))
)

save_model_plan(
  tcga_subtyped,
  "TCGA_pooled_clinical_adjusted",
  c("BUFFA_z", "SUBTYPE", tcga_pool_covars),
  diag_dir
)

rhs_tcga_pool_clin <- paste(c("BUFFA_z", "SUBTYPE", tcga_pool_covars), collapse = " + ")
m_tcga_pool_clinical <- run_cox_model(
  df = tcga_subtyped,
  rhs_formula = rhs_tcga_pool_clin,
  required_vars = c("BUFFA_z", "SUBTYPE", tcga_pool_covars),
  model_name = "TCGA_pooled_clinical_adjusted",
  out_dir = table_dir,
  relevel_list = list(STAGE_2 = "I_II", GRADE_2 = "G1_G2", ANY_TREATMENT = "No")
)

## ---------------------------------------------------------------------
## 5) TCGA subtype-specific analyses
## ---------------------------------------------------------------------

cat("\n=============================\nTCGA subtype-specific analyses\n=============================\n")

run_tcga_subtype_block <- function(df, subtype_name) {
  d <- df %>%
    filter(SUBTYPE == subtype_name) %>%
    mutate(hypoxia_group_st = factor(hypoxia_group_st, levels = c("High", "Low")))
  
  if (nrow(d) < 30) {
    save_skip_note(paste0("TCGA_", subtype_name, "_subtype_block"), "Subtype subset too small.", diag_dir)
    return(NULL)
  }
  
  covariate_availability(
    d,
    vars = c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"),
    file_stub = paste0("TCGA_", subtype_name)
  )
  
  file_stub <- paste0("TCGA_", subtype_name, "_withinSubtype")
  
  save_km_plot(
    df = d,
    formula_text = "Surv(time, event) ~ hypoxia_group_st",
    file_stub = file_stub,
    legend_title = paste0(subtype_name, " hypoxia")
  )
  
  save_model_plan(d, paste0(file_stub, "_cox_base"), c("hypoxia_group_st"), diag_dir)
  base_model <- run_cox_model(
    df = d,
    rhs_formula = "hypoxia_group_st",
    required_vars = c("hypoxia_group_st"),
    model_name = paste0(file_stub, "_cox_base"),
    out_dir = table_dir,
    relevel_list = list(hypoxia_group_st = "High")
  )
  
  subtype_main_covars <- intersect(
    c("AGE", "STAGE_2"),
    available_covars(d, c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"))
  )
  
  if (length(subtype_main_covars) > 0) {
    save_model_plan(
      d,
      paste0(file_stub, "_cox_clinical_main"),
      c("hypoxia_group_st", subtype_main_covars),
      diag_dir
    )
    
    clinical_main <- run_cox_model(
      df = d,
      rhs_formula = paste(c("hypoxia_group_st", subtype_main_covars), collapse = " + "),
      required_vars = c("hypoxia_group_st", subtype_main_covars),
      model_name = paste0(file_stub, "_cox_clinical_main"),
      out_dir = table_dir,
      relevel_list = list(hypoxia_group_st = "High", STAGE_2 = "I_II")
    )
  } else {
    clinical_main <- NULL
    save_skip_note(
      paste0(file_stub, "_cox_clinical_main"),
      "No usable AGE/STAGE_2 combination for subtype clinical model.",
      table_dir
    )
  }
  
  list(data = d, base = base_model, clinical_main = clinical_main)
}

tcga_luma  <- run_tcga_subtype_block(tcga_subtyped, "BRCA_LumA")
tcga_lumb  <- run_tcga_subtype_block(tcga_subtyped, "BRCA_LumB")
tcga_basal <- run_tcga_subtype_block(tcga_subtyped, "BRCA_Basal")

## ---------------------------------------------------------------------
## 6) TCGA Luminal B mutation enrichment
## ---------------------------------------------------------------------

cat("\n=============================\nTCGA LumB mutation enrichment\n=============================\n")

driver_genes <- c("TP53", "PIK3CA", "GATA3", "MAP3K1", "CDH1", "AKT1", "PTEN", "RB1", "NF1", "ERBB2", "FGFR1", "CCND1")

tcga_lumb_surv <- tcga_subtyped %>%
  filter(SUBTYPE == "BRCA_LumB") %>%
  distinct(SAMPLE_ID, PATIENT_ID, time, event, AGE, STAGE_2, GRADE_2, ANY_TREATMENT, hypoxia_group_st)

readr::write_csv(tcga_lumb_surv, file.path(data_dir, "TCGA_LumB_processed_survival_dataset.csv"))

mut_bin <- get_nonsyn_mut(tcga_mut, sample_ids = unique(tcga_lumb_surv$SAMPLE_ID))

if (!is.null(mut_bin)) {
  mut_drv <- mut_bin %>%
    filter(gene %in% driver_genes) %>%
    distinct(SAMPLE_ID, gene) %>%
    mutate(mut = 1)
  
  mat_drv <- mut_drv %>%
    tidyr::pivot_wider(names_from = gene, values_from = mut, values_fill = 0)
  
  dat_drv <- tcga_lumb_surv %>%
    select(SAMPLE_ID, hypoxia_group_st) %>%
    distinct() %>%
    left_join(mat_drv, by = "SAMPLE_ID")
  
  miss <- setdiff(driver_genes, names(dat_drv))
  if (length(miss) > 0) {
    for (g in miss) dat_drv[[g]] <- 0
  }
  
  dat_drv <- dat_drv %>%
    mutate(across(all_of(driver_genes), ~ tidyr::replace_na(., 0)))
  
  n_high <- sum(dat_drv$hypoxia_group_st == "High", na.rm = TRUE)
  n_low  <- sum(dat_drv$hypoxia_group_st == "Low", na.rm = TRUE)
  
  res_drv <- lapply(driver_genes, function(g) {
    tab <- table(dat_drv$hypoxia_group_st, factor(dat_drv[[g]], levels = c(0, 1)))
    colnames(tab) <- c("0", "1")
    out <- fisher_mutOR_high(tab)
    data.frame(
      gene = g,
      High_mut = tab["High", "1"], High_wt = tab["High", "0"],
      Low_mut  = tab["Low",  "1"], Low_wt  = tab["Low",  "0"],
      mutOR = out$mutOR, L95 = out$L95, U95 = out$U95, p = out$p
    )
  }) %>% bind_rows() %>%
    mutate(FDR = p.adjust(p, method = "BH")) %>%
    arrange(FDR, p)
  
  readr::write_csv(res_drv, file.path(table_dir, "TCGA_LumB_driver_mutation_enrichment.csv"))
  
  table4_out <- res_drv %>%
    mutate(
      High_mut_display = paste0(High_mut, "/", n_high),
      Low_mut_display  = paste0(Low_mut, "/", n_low),
      OR_round         = round(mutOR, 2),
      CI_display       = paste0(round(L95, 2), "–", round(U95, 2)),
      p                = signif(p, 3),
      FDR              = signif(FDR, 3)
    ) %>%
    select(gene, High_mut_display, Low_mut_display, OR_round, CI_display, p, FDR)
  
  readr::write_csv(table4_out, file.path(table_dir, "Table4_TP53_driver_enrichment.csv"))
}

## ---------------------------------------------------------------------
## 7) TCGA Luminal B CNA analyses
## ---------------------------------------------------------------------

cat("\n=============================\nTCGA LumB CNA analyses\n=============================\n")

tcga_cna_burden <- calc_cna_burden_continuous(tcga_cna)

if (!is.null(tcga_cna_burden)) {
  tcga_lumb_cna <- tcga_lumb_surv %>%
    left_join(tcga_cna_burden, by = "SAMPLE_ID") %>%
    filter(!is.na(CNA_mean_abs))
  
  if (IQR(tcga_lumb_cna$CNA_mean_abs, na.rm = TRUE) > 0) {
    tcga_lumb_cna <- tcga_lumb_cna %>%
      mutate(
        CNA_iqr2 = (CNA_mean_abs - median(CNA_mean_abs, na.rm = TRUE)) /
          IQR(CNA_mean_abs, na.rm = TRUE)
      )
  } else {
    tcga_lumb_cna <- tcga_lumb_cna %>% mutate(CNA_iqr2 = CNA_mean_abs)
  }
  
  readr::write_csv(tcga_lumb_cna, file.path(data_dir, "TCGA_LumB_CNA_processed_dataset.csv"))
  
  covariate_availability(
    tcga_lumb_cna,
    vars = c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT", "CNA_iqr2"),
    file_stub = "TCGA_LumB_CNA"
  )
  
  w_tcga_lumb_cna <- wilcox.test(CNA_mean_abs ~ hypoxia_group_st, data = tcga_lumb_cna)
  capture.output(w_tcga_lumb_cna, file = file.path(table_dir, "TCGA_LumB_CNA_wilcoxon.txt"))
  
  p_tcga_cna <- ggplot(tcga_lumb_cna, aes(x = hypoxia_group_st, y = CNA_mean_abs)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.35) +
    labs(
      title = paste0("TCGA LumB CNA burden by hypoxia (Wilcoxon p = ", signif(w_tcga_lumb_cna$p.value, 3), ")"),
      x = "Hypoxia group (within-subtype median)",
      y = "mean(|log2CNA|)"
    )
  ggsave(file.path(plot_dir, "TCGA_LumB_CNA_boxplot.png"), p_tcga_cna, width = 6, height = 4, dpi = 300)
  
  save_model_plan(tcga_lumb_cna, "TCGA_LumB_hypoxia_only", c("hypoxia_group_st"), diag_dir)
  m_tcga_lumb_hypoxia_only <- run_cox_model(
    df = tcga_lumb_cna,
    rhs_formula = "hypoxia_group_st",
    required_vars = c("hypoxia_group_st"),
    model_name = "TCGA_LumB_hypoxia_only",
    out_dir = table_dir,
    relevel_list = list(hypoxia_group_st = "High")
  )
  
  lumb_main_covars <- intersect(
    c("AGE", "STAGE_2"),
    available_covars(tcga_lumb_cna, c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT", "CNA_iqr2"))
  )
  
  save_model_plan(tcga_lumb_cna, "TCGA_LumB_hypoxia_clinical_main",
                  c("hypoxia_group_st", lumb_main_covars), diag_dir)
  m_tcga_lumb_clinical_main <- run_cox_model(
    df = tcga_lumb_cna,
    rhs_formula = paste(c("hypoxia_group_st", lumb_main_covars), collapse = " + "),
    required_vars = c("hypoxia_group_st", lumb_main_covars),
    model_name = "TCGA_LumB_hypoxia_clinical_main",
    out_dir = table_dir,
    relevel_list = list(hypoxia_group_st = "High", STAGE_2 = "I_II"),
    min_n = 50,
    min_events = 15
  )
  
  lumb_treat_covars <- intersect(c(lumb_main_covars, "ANY_TREATMENT"), names(tcga_lumb_cna))
  save_model_plan(tcga_lumb_cna, "TCGA_LumB_hypoxia_clinical_treatment_sensitivity",
                  c("hypoxia_group_st", lumb_treat_covars), diag_dir)
  m_tcga_lumb_clinical_treat <- run_cox_model(
    df = tcga_lumb_cna,
    rhs_formula = paste(c("hypoxia_group_st", lumb_treat_covars), collapse = " + "),
    required_vars = c("hypoxia_group_st", lumb_treat_covars),
    model_name = "TCGA_LumB_hypoxia_clinical_treatment_sensitivity",
    out_dir = table_dir,
    relevel_list = list(hypoxia_group_st = "High", STAGE_2 = "I_II", ANY_TREATMENT = "No"),
    min_n = 80,
    min_events = 15
  )
  
  save_model_plan(tcga_lumb_cna, "TCGA_LumB_hypoxia_clinical_CNA",
                  c("hypoxia_group_st", lumb_main_covars, "CNA_iqr2"), diag_dir)
  m_tcga_lumb_clinical_cna <- run_cox_model(
    df = tcga_lumb_cna,
    rhs_formula = paste(c("hypoxia_group_st", lumb_main_covars, "CNA_iqr2"), collapse = " + "),
    required_vars = c("hypoxia_group_st", lumb_main_covars, "CNA_iqr2"),
    model_name = "TCGA_LumB_hypoxia_clinical_CNA",
    out_dir = table_dir,
    relevel_list = list(hypoxia_group_st = "High", STAGE_2 = "I_II"),
    min_n = 50,
    min_events = 15
  )
  
  if (!is.null(m_tcga_lumb_clinical_main) && !is.null(m_tcga_lumb_clinical_cna)) {
    lrt_lumb <- anova(m_tcga_lumb_clinical_main$fit, m_tcga_lumb_clinical_cna$fit, test = "LRT")
    capture.output(lrt_lumb, file = file.path(table_dir, "TCGA_LumB_clinical_main_vs_clinicalCNA_LRT.txt"))
  }
}

## ---------------------------------------------------------------------
## 8) Load METABRIC
## ---------------------------------------------------------------------

cat("\n=============================\nLoading METABRIC\n=============================\n")

met_patient_path <- pick_first_file(metabric_dir, c("data_clinical_patient.txt"))
met_sample_path  <- pick_first_file(metabric_dir, c("data_clinical_sample.txt"))
met_expr_path <- pick_first_file(
  metabric_dir,
  c(
    "data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt",
    "data_mrna_illumina_microarray.txt",
    "data_mRNA.txt"
  )
)
met_cna_path     <- pick_first_file(metabric_dir, c("data_cna.txt"))
met_mut_path     <- pick_first_file(metabric_dir, c("data_mutations.txt"))

met_patient <- read_cbio_table(met_patient_path)
met_sample  <- read_cbio_table(met_sample_path)
met_expr    <- read_cbio_table(met_expr_path)
met_cna     <- if (!is.na(met_cna_path)) read_cbio_table(met_cna_path) else NULL
met_mut     <- if (!is.na(met_mut_path)) read_cbio_table(met_mut_path) else NULL

pid_col_pat <- if ("PATIENT_ID" %in% names(met_patient)) "PATIENT_ID" else names(met_patient)[1]
pid_col_sam <- if ("PATIENT_ID" %in% names(met_sample)) "PATIENT_ID" else NA_character_
sid_col_sam <- if ("SAMPLE_ID"  %in% names(met_sample)) "SAMPLE_ID"  else NA_character_
if (is.na(pid_col_sam) || is.na(sid_col_sam)) stop("METABRIC clinical_sample missing PATIENT_ID or SAMPLE_ID.")

met_master <- met_sample %>%
  mutate(
    PATIENT_ID = as.character(.data[[pid_col_sam]]),
    SAMPLE_ID  = as.character(.data[[sid_col_sam]])
  ) %>%
  left_join(
    met_patient %>% mutate(PATIENT_ID = as.character(.data[[pid_col_pat]])),
    by = "PATIENT_ID"
  )

hyp_genes <- c(
  "ALDOA", "ANGPTL4", "CA9", "ENO1", "HK2", "LDHA", "PGK1", "SLC2A1", "VEGFA",
  "PDK1", "ADM", "BNIP3", "NDRG1", "PFKFB3", "EGLN1", "EGLN3"
)

gene_col_expr <- pick_gene_col(met_expr)
if (is.na(gene_col_expr)) stop("Cannot identify gene symbol column in METABRIC expression table.")

drop_cols_expr <- unique(c(gene_col_expr, intersect(c("Entrez_Gene_Id"), names(met_expr))))
expr_sample_cols_all <- setdiff(names(met_expr), drop_cols_expr)
expr_sample_cols <- intersect(expr_sample_cols_all, met_master$SAMPLE_ID)

expr_sub <- met_expr %>%
  filter(.data[[gene_col_expr]] %in% hyp_genes) %>%
  select(all_of(gene_col_expr), all_of(expr_sample_cols))

expr_sub2 <- expr_sub %>%
  mutate(GENE = as.character(.data[[gene_col_expr]])) %>%
  group_by(GENE) %>%
  summarise(across(all_of(expr_sample_cols), ~ mean(to_num(.), na.rm = TRUE)), .groups = "drop")

mat_expr <- as.matrix(expr_sub2[, expr_sample_cols, drop = FALSE])
storage.mode(mat_expr) <- "numeric"

z_expr <- t(scale(t(mat_expr)))
z_expr[!is.finite(z_expr)] <- NA
hyp_score <- colMeans(z_expr, na.rm = TRUE)

met_hyp_df <- data.frame(
  SAMPLE_ID = names(hyp_score),
  HYPOXIA_SCORE = as.numeric(hyp_score)
)

met_master <- met_master %>%
  mutate(SAMPLE_ID = as.character(SAMPLE_ID)) %>%
  left_join(met_hyp_df, by = "SAMPLE_ID")

met_surv_cols   <- pick_survival_cols(met_master)
met_subtype_col <- pick_subtype_col_metabric(met_master)
met_age_col     <- pick_age_col(met_master)
met_stage_col   <- pick_stage_col(met_master)
met_grade_col   <- pick_grade_col(met_master)

met_df <- met_master %>%
  transmute(
    PATIENT_ID    = as.character(PATIENT_ID),
    SAMPLE_ID     = as.character(SAMPLE_ID),
    time          = to_num(.data[[met_surv_cols$time]]),
    event         = status_to_event(.data[[met_surv_cols$status]]),
    HYPOXIA_SCORE = to_num(HYPOXIA_SCORE),
    HYPOXIA_Z     = as.numeric(scale(HYPOXIA_SCORE)),
    SUBTYPE_ANY   = if (!is.na(met_subtype_col)) as.character(.data[[met_subtype_col]]) else NA_character_,
    AGE           = if (!is.na(met_age_col)) to_num(.data[[met_age_col]]) else NA_real_,
    STAGE_RAW     = if (!is.na(met_stage_col)) as.character(.data[[met_stage_col]]) else NA_character_,
    GRADE_RAW     = if (!is.na(met_grade_col)) as.character(.data[[met_grade_col]]) else NA_character_
  ) %>%
  mutate(
    STAGE_SIMPLE = simplify_stage(STAGE_RAW),
    STAGE_2      = stage_to_2cat(STAGE_SIMPLE),
    GRADE_SIMPLE = simplify_grade(GRADE_RAW),
    GRADE_2      = grade_to_2cat(GRADE_SIMPLE)
  )

met_df$ANY_TREATMENT <- derive_any_treatment(met_master)

met_df <- met_df %>%
  filter(!is.na(time), time > 0, !is.na(event), !is.na(HYPOXIA_SCORE)) %>%
  mutate(
    hypoxia_group = ifelse(HYPOXIA_SCORE >= median(HYPOXIA_SCORE, na.rm = TRUE), "High", "Low"),
    hypoxia_group = factor(hypoxia_group, levels = c("High", "Low"))
  )

readr::write_csv(met_df, file.path(data_dir, "METABRIC_processed_survival_dataset.csv"))

cat("METABRIC rows:", nrow(met_df), "\n")
cat("METABRIC events:", sum(met_df$event == 1, na.rm = TRUE), "\n")
cat("METABRIC endpoint:", met_surv_cols$endpoint, "\n")

covariate_availability(
  met_df,
  vars = c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"),
  file_stub = "METABRIC_overall"
)

## ---------------------------------------------------------------------
## 9) METABRIC survival analyses
## ---------------------------------------------------------------------

cat("\n=============================\nMETABRIC survival analyses\n=============================\n")

save_km_plot(
  df = met_df,
  formula_text = "Surv(time, event) ~ hypoxia_group",
  file_stub = "METABRIC_KM_global_median",
  legend_title = "METABRIC hypoxia"
)

save_model_plan(met_df, "METABRIC_hypoxia_only", c("hypoxia_group"), diag_dir)
m_met_base <- run_cox_model(
  df = met_df,
  rhs_formula = "hypoxia_group",
  required_vars = c("hypoxia_group"),
  model_name = "METABRIC_hypoxia_only",
  out_dir = table_dir,
  relevel_list = list(hypoxia_group = "High")
)

met_subtyped <- met_df %>%
  filter(!is.na(SUBTYPE_ANY), SUBTYPE_ANY != "") %>%
  mutate(SUBTYPE_ANY = factor(SUBTYPE_ANY))

save_model_plan(met_subtyped, "METABRIC_subtype_adjusted", c("hypoxia_group", "SUBTYPE_ANY"), diag_dir)
m_met_subtype <- run_cox_model(
  df = met_subtyped,
  rhs_formula = "hypoxia_group + SUBTYPE_ANY",
  required_vars = c("hypoxia_group", "SUBTYPE_ANY"),
  model_name = "METABRIC_subtype_adjusted",
  out_dir = table_dir,
  relevel_list = list(hypoxia_group = "High")
)

met_pool_avail <- covariate_availability(
  met_subtyped,
  vars = c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"),
  file_stub = "METABRIC_pooled"
)

met_pool_covars <- intersect(
  c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"),
  available_covars(met_subtyped, c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"))
)

save_model_plan(met_subtyped, "METABRIC_clinical_adjusted",
                c("hypoxia_group", "SUBTYPE_ANY", met_pool_covars), diag_dir)
m_met_clinical <- run_cox_model(
  df = met_subtyped,
  rhs_formula = paste(c("hypoxia_group", "SUBTYPE_ANY", met_pool_covars), collapse = " + "),
  required_vars = c("hypoxia_group", "SUBTYPE_ANY", met_pool_covars),
  model_name = "METABRIC_clinical_adjusted",
  out_dir = table_dir,
  relevel_list = list(hypoxia_group = "High", STAGE_2 = "I_II", GRADE_2 = "G1_G2", ANY_TREATMENT = "No")
)

save_model_plan(met_subtyped, "METABRIC_clinical_adjusted_continuousHypoxia",
                c("HYPOXIA_Z", "SUBTYPE_ANY", met_pool_covars), diag_dir)
m_met_clinical_cont <- run_cox_model(
  df = met_subtyped,
  rhs_formula = paste(c("HYPOXIA_Z", "SUBTYPE_ANY", met_pool_covars), collapse = " + "),
  required_vars = c("HYPOXIA_Z", "SUBTYPE_ANY", met_pool_covars),
  model_name = "METABRIC_clinical_adjusted_continuousHypoxia",
  out_dir = table_dir,
  relevel_list = list(STAGE_2 = "I_II", GRADE_2 = "G1_G2", ANY_TREATMENT = "No")
)

## ---------------------------------------------------------------------
## 10) METABRIC CNA analyses
## ---------------------------------------------------------------------

cat("\n=============================\nMETABRIC CNA analyses\n=============================\n")

met_cna_burden <- calc_cna_burden_discrete(met_cna)

if (!is.null(met_cna_burden)) {
  met_cna_df <- met_df %>%
    left_join(met_cna_burden, by = "SAMPLE_ID") %>%
    filter(!is.na(CNA_mean_abs))
  
  readr::write_csv(met_cna_df, file.path(data_dir, "METABRIC_CNA_processed_dataset.csv"))
  
  w_met_cna <- wilcox.test(CNA_mean_abs ~ hypoxia_group, data = met_cna_df)
  capture.output(w_met_cna, file = file.path(table_dir, "METABRIC_CNA_wilcoxon.txt"))
  
  p_met_cna <- ggplot(met_cna_df, aes(x = hypoxia_group, y = CNA_mean_abs)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.30) +
    labs(
      title = paste0("METABRIC CNA burden by hypoxia (Wilcoxon p = ", signif(w_met_cna$p.value, 3), ")"),
      x = "Hypoxia group",
      y = "mean(|discrete CNA|)"
    )
  ggsave(file.path(plot_dir, "METABRIC_CNA_boxplot.png"), p_met_cna, width = 6, height = 4, dpi = 300)
  
  met_cna_subtyped <- met_cna_df %>%
    filter(!is.na(SUBTYPE_ANY), SUBTYPE_ANY != "") %>%
    mutate(SUBTYPE_ANY = factor(SUBTYPE_ANY))
  
  covariate_availability(
    met_cna_subtyped,
    vars = c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT", "CNA_mean_abs"),
    file_stub = "METABRIC_CNA"
  )
  
  met_cna_covars <- intersect(
    c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"),
    available_covars(met_cna_subtyped, c("AGE", "STAGE_2", "GRADE_2", "ANY_TREATMENT"))
  )
  
  save_model_plan(met_cna_subtyped, "METABRIC_subtypeClinical_hypoxia",
                  c("hypoxia_group", "SUBTYPE_ANY", met_cna_covars), diag_dir)
  m_met_cna_base <- run_cox_model(
    df = met_cna_subtyped,
    rhs_formula = paste(c("hypoxia_group", "SUBTYPE_ANY", met_cna_covars), collapse = " + "),
    required_vars = c("hypoxia_group", "SUBTYPE_ANY", met_cna_covars),
    model_name = "METABRIC_subtypeClinical_hypoxia",
    out_dir = table_dir,
    relevel_list = list(hypoxia_group = "High", STAGE_2 = "I_II", GRADE_2 = "G1_G2", ANY_TREATMENT = "No")
  )
  
  save_model_plan(met_cna_subtyped, "METABRIC_subtypeClinical_hypoxia_CNA",
                  c("hypoxia_group", "SUBTYPE_ANY", met_cna_covars, "CNA_mean_abs"), diag_dir)
  m_met_cna_full <- run_cox_model(
    df = met_cna_subtyped,
    rhs_formula = paste(c("hypoxia_group", "SUBTYPE_ANY", met_cna_covars, "CNA_mean_abs"), collapse = " + "),
    required_vars = c("hypoxia_group", "SUBTYPE_ANY", met_cna_covars, "CNA_mean_abs"),
    model_name = "METABRIC_subtypeClinical_hypoxia_CNA",
    out_dir = table_dir,
    relevel_list = list(hypoxia_group = "High", STAGE_2 = "I_II", GRADE_2 = "G1_G2", ANY_TREATMENT = "No")
  )
  
  if (!is.null(m_met_cna_base) && !is.null(m_met_cna_full)) {
    lrt_met <- anova(m_met_cna_base$fit, m_met_cna_full$fit, test = "LRT")
    capture.output(lrt_met, file = file.path(table_dir, "METABRIC_subtypeClinical_hypoxia_vs_hypoxiaCNA_LRT.txt"))
  }
}

## ---------------------------------------------------------------------
## 11) Model summary tables
## ---------------------------------------------------------------------

cat("\n=============================\nWriting model summary tables\n=============================\n")

analysis_summary <- list()

analysis_summary[[1]] <- data.frame(
  analysis = "TCGA pooled clinical-adjusted",
  covariates_used = paste(c("BUFFA_z", "SUBTYPE", tcga_pool_covars), collapse = " + ")
)

analysis_summary[[2]] <- data.frame(
  analysis = "METABRIC pooled clinical-adjusted",
  covariates_used = paste(c("hypoxia_group", "SUBTYPE_ANY", met_pool_covars), collapse = " + ")
)

lumb_main_string <- if (exists("lumb_main_covars")) {
  paste(c("hypoxia_group_st", lumb_main_covars), collapse = " + ")
} else {
  "Not run"
}

analysis_summary[[3]] <- data.frame(
  analysis = "TCGA LumB clinical main",
  covariates_used = lumb_main_string
)

analysis_summary_tab <- bind_rows(analysis_summary)
readr::write_csv(analysis_summary_tab, file.path(diag_dir, "analysis_model_summary.csv"))

## ---------------------------------------------------------------------
## 12) Reproducibility manifest
## ---------------------------------------------------------------------

cat("\n=============================\nWriting reproducibility manifest\n=============================\n")

manifest <- data.frame(
  cohort = c("TCGA", "METABRIC"),
  dataset_folder = c(tcga_dir, metabric_dir),
  survival_endpoint = c(tcga_surv_cols$endpoint, met_surv_cols$endpoint),
  processed_file = c(
    "TCGA_processed_survival_dataset.csv",
    "METABRIC_processed_survival_dataset.csv"
  )
)
readr::write_csv(manifest, file.path(out_dir, "dataset_manifest.csv"))

if (length(model_registry) > 0) {
  readr::write_csv(bind_rows(model_registry), file.path(diag_dir, "all_model_registry.csv"))
}
if (length(ph_registry) > 0) {
  readr::write_csv(bind_rows(ph_registry), file.path(diag_dir, "all_PH_tests_combined.csv"))
}
if (length(skip_registry) > 0) {
  readr::write_csv(bind_rows(skip_registry), file.path(diag_dir, "all_skipped_models.csv"))
}

session_txt <- capture.output(sessionInfo())
writeLines(session_txt, file.path(out_dir, "sessionInfo.txt"))

notes <- c(
  "This script implements multivariable adjustment using standard clinical covariates.",
  "TCGA pooled model uses BUFFA_z + SUBTYPE + all usable standard clinical covariates.",
  "METABRIC pooled model uses hypoxia + SUBTYPE + all usable standard clinical covariates.",
  "METABRIC stage/grade recoding handles numeric coding, so STAGE_2 and GRADE_2 are retained when available.",
  "Subtype-specific TCGA models use parsimonious AGE + STAGE_2 adjustment because sparse event counts make overfitting a real risk.",
  "TCGA Luminal B includes separate sensitivity / extension models for treatment and CNA burden.",
  "Every Cox model exports PH diagnostics as *_PH_test.csv and *_PH_plots.pdf.",
  "Processed per-sample analysis files are exported under outputs_final/processed_data/.",
  "Important: use the newly generated output files as the source of final manuscript numbers."
)
writeLines(notes, file.path(out_dir, "README_analysis_pipeline.txt"))

cat("\nAll primary outputs are in:\n", out_dir, "\n")

## ---------------------------------------------------------------------
## 13) Final model output summary
## ---------------------------------------------------------------------

cat("\n=============================\nCreating final model outputs\n=============================\n")

read_csv_safe <- function(path) {
  if (!file.exists(path)) return(NULL)
  suppressMessages(readr::read_csv(path, show_col_types = FALSE))
}

model_info_df <- if (length(model_registry) > 0) bind_rows(model_registry) else NULL
ph_df_all     <- if (length(ph_registry) > 0) bind_rows(ph_registry) else NULL

tidy_files <- list.files(table_dir, pattern = "_tidy\\.csv$", full.names = TRUE)
model_terms_df <- if (length(tidy_files) > 0) {
  purrr::map_dfr(tidy_files, function(fp) {
    x <- read_csv_safe(fp)
    if (is.null(x)) return(NULL)
    x$model <- sub("_tidy\\.csv$", "", basename(fp))
    x
  })
} else NULL

cox_summary_df <- if (!is.null(model_terms_df) && !is.null(model_info_df)) {
  model_terms_df %>%
    left_join(model_info_df, by = "model") %>%
    select(model, formula, n, events, term, HR, CI_low, CI_high, p.value, everything())
} else {
  model_terms_df
}

complete_case_df <- {
  files <- list.files(diag_dir, pattern = "_planned_complete_cases\\.csv$", full.names = TRUE)
  if (length(files) > 0) {
    purrr::map_dfr(files, function(fp) {
      x <- read_csv_safe(fp)
      if (is.null(x)) return(NULL)
      x$analysis <- sub("_planned_complete_cases\\.csv$", "", basename(fp))
      x
    })
  } else NULL
}

covariate_avail_df <- {
  files <- list.files(diag_dir, pattern = "_covariate_availability\\.csv$", full.names = TRUE)
  if (length(files) > 0) {
    purrr::map_dfr(files, function(fp) {
      x <- read_csv_safe(fp)
      if (is.null(x)) return(NULL)
      x$analysis <- sub("_covariate_availability\\.csv$", "", basename(fp))
      x
    })
  } else NULL
}

analysis_model_df <- read_csv_safe(file.path(diag_dir, "analysis_model_summary.csv"))
mut_df            <- read_csv_safe(file.path(table_dir, "TCGA_LumB_driver_mutation_enrichment.csv"))
manifest_df       <- read_csv_safe(file.path(out_dir, "dataset_manifest.csv"))

fit_weibull_aft <- function(df, rhs_formula, required_vars, model_name,
                            relevel_list = list(), min_n = 30, min_events = 10) {
  keep_vars <- unique(c("time", "event", required_vars))
  keep_vars <- keep_vars[keep_vars %in% names(df)]
  
  d <- df %>%
    dplyr::select(all_of(keep_vars)) %>%
    filter(complete.cases(.), time > 0)
  
  if (nrow(d) < min_n) return(NULL)
  if (sum(d$event == 1, na.rm = TRUE) < min_events) return(NULL)
  
  for (nm in names(relevel_list)) {
    if (nm %in% names(d)) {
      d[[nm]] <- factor(d[[nm]])
      ref <- relevel_list[[nm]]
      if (ref %in% levels(d[[nm]])) {
        d[[nm]] <- relevel(d[[nm]], ref = ref)
      }
    }
  }
  
  ftxt <- paste0("Surv(time, event) ~ ", rhs_formula)
  fit <- survival::survreg(as.formula(ftxt), data = d, dist = "weibull")
  
  sm <- summary(fit)
  coef_tab <- as.data.frame(sm$table)
  coef_tab$term <- rownames(coef_tab)
  rownames(coef_tab) <- NULL
  
  names(coef_tab) <- sub("Value", "estimate", names(coef_tab))
  names(coef_tab) <- sub("Std\\. Error", "std.error", names(coef_tab))
  names(coef_tab) <- sub("z", "statistic", names(coef_tab))
  names(coef_tab) <- sub("p", "p.value", names(coef_tab))
  
  ci <- suppressMessages(confint(fit))
  ci_df <- data.frame(term = rownames(ci), conf.low = ci[, 1], conf.high = ci[, 2], row.names = NULL)
  
  out <- coef_tab %>%
    left_join(ci_df, by = "term") %>%
    mutate(
      TimeRatio = exp(estimate),
      TR_low    = exp(conf.low),
      TR_high   = exp(conf.high),
      model     = model_name,
      formula   = ftxt,
      n         = nrow(d),
      events    = sum(d$event == 1, na.rm = TRUE)
    )
  
  list(fit = fit, tidy = out, data = d)
}

pick_term_row <- function(df, model_name, term_pattern) {
  if (is.null(df)) return(NULL)
  x <- df %>% filter(model == model_name, grepl(term_pattern, term, ignore.case = TRUE))
  if (nrow(x) == 0) return(NULL)
  x[1, , drop = FALSE]
}

get_ph_p <- function(ph_df, model_name, term_exact) {
  if (is.null(ph_df)) return(NA_real_)
  x <- ph_df %>% filter(model == model_name, term == term_exact)
  if (nrow(x) == 0) return(NA_real_)
  as.numeric(x$p[1])
}

fmt_num <- function(x, k = 3) {
  ifelse(is.na(x), "NA", format(round(x, k), nsmall = k, trim = TRUE))
}

fmt_p <- function(x) {
  ifelse(is.na(x), "NA", format(signif(x, 3), scientific = TRUE))
}

aft_specs <- list(
  list(
    model = "METABRIC_hypoxia_only_AFT",
    data_obj = met_df,
    rhs = "hypoxia_group",
    vars = c("hypoxia_group"),
    relevel = list(hypoxia_group = "High")
  ),
  list(
    model = "METABRIC_clinical_adjusted_AFT",
    data_obj = met_subtyped,
    rhs = paste(c("hypoxia_group", "SUBTYPE_ANY", met_pool_covars), collapse = " + "),
    vars = c("hypoxia_group", "SUBTYPE_ANY", met_pool_covars),
    relevel = list(hypoxia_group = "High", STAGE_2 = "I_II", GRADE_2 = "G1_G2", ANY_TREATMENT = "No")
  )
)

if (exists("met_cna_subtyped") && exists("met_cna_covars")) {
  aft_specs[[length(aft_specs) + 1L]] <- list(
    model = "METABRIC_subtypeClinical_hypoxia_CNA_AFT",
    data_obj = met_cna_subtyped,
    rhs = paste(c("hypoxia_group", "SUBTYPE_ANY", met_cna_covars, "CNA_mean_abs"), collapse = " + "),
    vars = c("hypoxia_group", "SUBTYPE_ANY", met_cna_covars, "CNA_mean_abs"),
    relevel = list(hypoxia_group = "High", STAGE_2 = "I_II", GRADE_2 = "G1_G2", ANY_TREATMENT = "No")
  )
}

aft_results <- purrr::map(aft_specs, function(sp) {
  fit_weibull_aft(
    df = sp$data_obj,
    rhs_formula = sp$rhs,
    required_vars = sp$vars,
    model_name = sp$model,
    relevel_list = sp$relevel
  )
})

aft_summary_df <- purrr::map_dfr(aft_results, function(x) {
  if (is.null(x)) return(NULL)
  x$tidy
})

build_preferred_row <- function(label, cox_model, cox_term_pattern, ph_term_exact, aft_model = NULL) {
  cox_row <- pick_term_row(cox_summary_df, cox_model, cox_term_pattern)
  
  if (is.null(cox_row)) {
    return(data.frame(
      label = label,
      preferred_method = "Not available",
      effect_type = NA_character_,
      estimate = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p.value = NA_real_,
      n = NA_real_,
      events = NA_real_,
      cox_model = cox_model,
      cox_ph_global_p = NA_real_,
      cox_ph_exposure_p = NA_real_,
      note = "Result not available"
    ))
  }
  
  ph_global_p   <- get_ph_p(ph_df_all, cox_model, "GLOBAL")
  ph_exposure_p <- get_ph_p(ph_df_all, cox_model, ph_term_exact)
  ph_ok <- !is.na(ph_global_p) && !is.na(ph_exposure_p) && ph_global_p >= 0.05 && ph_exposure_p >= 0.05
  
  if (ph_ok || is.null(aft_model)) {
    return(data.frame(
      label = label,
      preferred_method = "Cox PH model",
      effect_type = "Hazard ratio",
      estimate = cox_row$HR,
      ci_low = cox_row$CI_low,
      ci_high = cox_row$CI_high,
      p.value = cox_row$p.value,
      n = cox_row$n,
      events = cox_row$events,
      cox_model = cox_model,
      cox_ph_global_p = ph_global_p,
      cox_ph_exposure_p = ph_exposure_p,
      note = "Cox HR is reported with model-specific PH diagnostics."
    ))
  }
  
  aft_row <- pick_term_row(aft_summary_df, aft_model, cox_term_pattern)
  if (is.null(aft_row)) {
    return(data.frame(
      label = label,
      preferred_method = "Cox PH model",
      effect_type = "Hazard ratio",
      estimate = cox_row$HR,
      ci_low = cox_row$CI_low,
      ci_high = cox_row$CI_high,
      p.value = cox_row$p.value,
      n = cox_row$n,
      events = cox_row$events,
      cox_model = cox_model,
      cox_ph_global_p = ph_global_p,
      cox_ph_exposure_p = ph_exposure_p,
      note = "PH diagnostics suggested violation, but the Weibull AFT model was not available; Cox HR is retained with caution."
    ))
  }
  
  data.frame(
    label = label,
    preferred_method = "Weibull AFT model",
    effect_type = "Time ratio",
    estimate = aft_row$TimeRatio,
    ci_low = aft_row$TR_low,
    ci_high = aft_row$TR_high,
    p.value = aft_row$p.value,
    n = aft_row$n,
    events = aft_row$events,
    cox_model = cox_model,
    cox_ph_global_p = ph_global_p,
    cox_ph_exposure_p = ph_exposure_p,
    note = "Cox PH violated; prefer Weibull AFT time ratio for interpretation."
  )
}

preferred_results_df <- bind_rows(
  build_preferred_row("TCGA pooled subtype-adjusted hypoxia", "TCGA_pooled_subtype_adjusted", "^BUFFA_z$", "BUFFA_z"),
  build_preferred_row("TCGA pooled clinical-adjusted hypoxia", "TCGA_pooled_clinical_adjusted", "^BUFFA_z$", "BUFFA_z"),
  build_preferred_row("TCGA LumB base hypoxia", "TCGA_BRCA_LumB_withinSubtype_cox_base", "^hypoxia_group_st", "hypoxia_group_st"),
  build_preferred_row("TCGA LumB clinical-main hypoxia", "TCGA_BRCA_LumB_withinSubtype_cox_clinical_main", "^hypoxia_group_st", "hypoxia_group_st"),
  build_preferred_row("TCGA LumB clinical+CNA hypoxia", "TCGA_LumB_hypoxia_clinical_CNA", "^hypoxia_group_st", "hypoxia_group_st"),
  build_preferred_row("TCGA LumB clinical+CNA burden", "TCGA_LumB_hypoxia_clinical_CNA", "^CNA", "CNA_iqr2"),
  build_preferred_row("METABRIC hypoxia-only", "METABRIC_hypoxia_only", "^hypoxia_group", "hypoxia_group", "METABRIC_hypoxia_only_AFT"),
  build_preferred_row("METABRIC clinical-adjusted hypoxia", "METABRIC_clinical_adjusted", "^hypoxia_group", "hypoxia_group", "METABRIC_clinical_adjusted_AFT"),
  build_preferred_row("METABRIC subtype+clinical+CNA hypoxia", "METABRIC_subtypeClinical_hypoxia_CNA", "^hypoxia_group", "hypoxia_group", "METABRIC_subtypeClinical_hypoxia_CNA_AFT"),
  build_preferred_row("METABRIC subtype+clinical+CNA burden", "METABRIC_subtypeClinical_hypoxia_CNA", "^CNA_mean_abs$", "CNA_mean_abs", "METABRIC_subtypeClinical_hypoxia_CNA_AFT")
)

ph_flag_df <- if (!is.null(ph_df_all)) {
  key_map <- tibble::tribble(
    ~label, ~model, ~exposure_term,
    "TCGA pooled subtype-adjusted hypoxia", "TCGA_pooled_subtype_adjusted", "BUFFA_z",
    "TCGA pooled clinical-adjusted hypoxia", "TCGA_pooled_clinical_adjusted", "BUFFA_z",
    "TCGA LumB base hypoxia", "TCGA_BRCA_LumB_withinSubtype_cox_base", "hypoxia_group_st",
    "TCGA LumB clinical-main hypoxia", "TCGA_BRCA_LumB_withinSubtype_cox_clinical_main", "hypoxia_group_st",
    "TCGA LumB clinical+CNA hypoxia", "TCGA_LumB_hypoxia_clinical_CNA", "hypoxia_group_st",
    "TCGA LumB clinical+CNA burden", "TCGA_LumB_hypoxia_clinical_CNA", "CNA_iqr2",
    "METABRIC hypoxia-only", "METABRIC_hypoxia_only", "hypoxia_group",
    "METABRIC clinical-adjusted hypoxia", "METABRIC_clinical_adjusted", "hypoxia_group",
    "METABRIC subtype+clinical+CNA hypoxia", "METABRIC_subtypeClinical_hypoxia_CNA", "hypoxia_group",
    "METABRIC subtype+clinical+CNA burden", "METABRIC_subtypeClinical_hypoxia_CNA", "CNA_mean_abs"
  )
  
  key_map %>%
    rowwise() %>%
    mutate(
      global_p = get_ph_p(ph_df_all, model, "GLOBAL"),
      exposure_p = get_ph_p(ph_df_all, model, exposure_term),
      ph_ok = !is.na(global_p) & !is.na(exposure_p) & global_p >= 0.05 & exposure_p >= 0.05
    ) %>%
    ungroup()
} else NULL

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "preferred_results")
openxlsx::writeData(wb, "preferred_results", preferred_results_df)

openxlsx::addWorksheet(wb, "cox_results")
openxlsx::writeData(wb, "cox_results", if (is.null(cox_summary_df)) data.frame(message = "No Cox model summary found") else cox_summary_df)

openxlsx::addWorksheet(wb, "PH_summary")
openxlsx::writeData(wb, "PH_summary", if (is.null(ph_flag_df)) data.frame(message = "No PH diagnostics found") else ph_flag_df)

openxlsx::addWorksheet(wb, "AFT_models")
openxlsx::writeData(wb, "AFT_models", if (is.null(aft_summary_df) || nrow(aft_summary_df) == 0) data.frame(message = "No Weibull AFT models were run") else aft_summary_df)

openxlsx::addWorksheet(wb, "complete_cases")
openxlsx::writeData(wb, "complete_cases", if (is.null(complete_case_df)) data.frame(message = "No complete-case table found") else complete_case_df)

openxlsx::addWorksheet(wb, "covariate_availability")
openxlsx::writeData(wb, "covariate_availability", if (is.null(covariate_avail_df)) data.frame(message = "No covariate availability table found") else covariate_avail_df)

openxlsx::addWorksheet(wb, "analysis_model_plan")
openxlsx::writeData(wb, "analysis_model_plan", if (is.null(analysis_model_df)) data.frame(message = "No model summary found") else analysis_model_df)

openxlsx::addWorksheet(wb, "TCGA_LumB_mutation")
openxlsx::writeData(wb, "TCGA_LumB_mutation", if (is.null(mut_df)) data.frame(message = "No mutation enrichment output found") else mut_df)

openxlsx::addWorksheet(wb, "manifest")
openxlsx::writeData(wb, "manifest", if (is.null(manifest_df)) data.frame(message = "No manifest found") else manifest_df)

openxlsx::saveWorkbook(wb, file.path(out_dir, "00_final_results_PHsafe.xlsx"), overwrite = TRUE)

make_summary_line <- function(label_pattern) {
  x <- preferred_results_df %>% filter(label == label_pattern)
  if (nrow(x) == 0) return(paste0(label_pattern, ": not available"))
  x <- x[1, ]
  
  paste0(
    x$label, ": ",
    x$preferred_method, " -> ",
    ifelse(x$effect_type == "Hazard ratio", "HR=", "TR="),
    fmt_num(x$estimate),
    " (", fmt_num(x$ci_low), "-", fmt_num(x$ci_high), ")",
    ", p=", fmt_p(x$p.value),
    ", n=", x$n, ", events=", x$events,
    ", Cox PH global p=", fmt_p(x$cox_ph_global_p),
    ", exposure PH p=", fmt_p(x$cox_ph_exposure_p)
  )
}

summary_lines <- c(
  "PH-safe final output summary",
  "Output files are provided in this repository.",
  "",
  make_summary_line("TCGA pooled subtype-adjusted hypoxia"),
  make_summary_line("TCGA pooled clinical-adjusted hypoxia"),
  make_summary_line("TCGA LumB base hypoxia"),
  make_summary_line("TCGA LumB clinical-main hypoxia"),
  make_summary_line("TCGA LumB clinical+CNA hypoxia"),
  make_summary_line("TCGA LumB clinical+CNA burden"),
  make_summary_line("METABRIC hypoxia-only"),
  make_summary_line("METABRIC clinical-adjusted hypoxia"),
  make_summary_line("METABRIC subtype+clinical+CNA hypoxia"),
  make_summary_line("METABRIC subtype+clinical+CNA burden"),
  "",
  "Open 00_final_results_PHsafe.xlsx first.",
  "Sheet preferred_results contains the main numbers to cite.",
  "For METABRIC models with proportional hazards violations, Weibull AFT time ratios are used as the preferred interpretable estimates."
)
writeLines(summary_lines, file.path(out_dir, "00_final_summary_PHsafe.txt"))

writeLines(c(
  "This output folder contains the final files used for reporting and reproducibility.",
  "Open first: 00_final_results_PHsafe.xlsx",
  "Then read: 00_final_summary_PHsafe.txt",
  "Key processed datasets are kept under processed_data/",
  "Main figures are kept under plots/",
  "Detailed per-model tables and diagnostic files are retained under tables/ and diagnostics/."
), file.path(out_dir, "README_outputs.txt"))

key_plot_files <- c(
  file.path(plot_dir, "TCGA_KM_global_median.png"),
  file.path(plot_dir, "TCGA_global_hypoxia_by_subtype_barplot.png"),
  file.path(plot_dir, "TCGA_BRCA_LumB_withinSubtype.png"),
  file.path(plot_dir, "TCGA_LumB_CNA_boxplot.png"),
  file.path(plot_dir, "METABRIC_KM_global_median.png"),
  file.path(plot_dir, "METABRIC_CNA_boxplot.png")
)

keep_files <- c(
  file.path(out_dir, "00_final_results_PHsafe.xlsx"),
  file.path(out_dir, "00_final_summary_PHsafe.txt"),
  file.path(out_dir, "README_outputs.txt"),
  file.path(out_dir, "dataset_manifest.csv"),
  file.path(out_dir, "sessionInfo.txt"),
  file.path(data_dir, "TCGA_processed_survival_dataset.csv"),
  file.path(data_dir, "METABRIC_processed_survival_dataset.csv"),
  file.path(data_dir, "TCGA_LumB_CNA_processed_dataset.csv"),
  file.path(data_dir, "METABRIC_CNA_processed_dataset.csv"),
  key_plot_files
)

keep_files <- keep_files[file.exists(keep_files)]
all_files  <- list.files(out_dir, recursive = TRUE, full.names = TRUE)

## Intermediate tables and diagnostic files are retained so that the model
## outputs can be inspected if needed. The objects in keep_files are the core
## files used for manuscript reporting.

cat("\nFinal outputs written to:\n", out_dir, "\n")
cat("Open this first:\n", file.path(out_dir, "00_final_results_PHsafe.xlsx"), "\n")

library(dplyr)
library(tidyr)
library(stringr)

# --- 1. Read the export ------------------------------------------------------
df <- read.csv("~/Downloads/NIHTB_corrected_raw.csv", stringsAsFactors = FALSE)

score_cols <- c("nsm_rawscore", "arw_ratescore", "scs_ratescore")
# mfs_rawscore (Sequences) and arp_orderedpairs (Arranging Pictures) dropped
# per PI direction.

date_cols <- c("startdate", "startdate_2", "startdate_3", "startdate_4", "startdate_5")

coalesce_first <- function(x) {
  non_na <- x[!is.na(x) & x != ""]
  if (length(non_na) == 0) NA else non_na[1]
}

# --- 2. Sanity check BEFORE collapsing ---------------------------------------
conflicts <- df %>%
  group_by(record_id, redcap_event_name) %>%
  summarise(across(all_of(score_cols), ~ sum(!is.na(.))), .groups = "drop") %>%
  filter(if_any(all_of(score_cols), ~ . > 1))

if (nrow(conflicts) > 0) {
  warning("Duplicate administrations found - review `conflicts` before trusting the output.")
  print(conflicts)
}

# --- 3. Collapse to one row per record x visit -------------------------------
clean <- df %>%
  filter(!is.na(redcap_repeat_instrument) & redcap_repeat_instrument != "") %>%
  group_by(record_id, redcap_event_name) %>%
  summarise(across(all_of(score_cols), coalesce_first), .groups = "drop") %>%
  mutate(
    visit_num = as.numeric(str_extract(redcap_event_name, "(?<=v)\\d+")),
    visit = paste0("v", visit_num)
  ) %>%
  select(record_id, visit, all_of(score_cols))

# --- 4. Protocol % change composite, v1 (baseline) -> v6 --------------------
# % change per subtest = (follow-up - baseline) / baseline * 100
# Composite = mean of the 5 subtest % changes.
# Missing data rule per protocol: composite requires >= 4 of 5 valid subtest
# scores at BOTH visits (a subtest can't contribute a % change if it's missing
# at either baseline or follow-up). Fewer than 4 -> "unevaluable", no AE call.
wide <- clean %>%
  filter(visit %in% c("v1", "v6")) %>%
  pivot_wider(id_cols = record_id, names_from = visit, values_from = all_of(score_cols))

# Flag zero/near-zero baselines - % change is undefined (division by zero) or
# explosively unstable for tiny denominators. Check this before trusting any
# % change values below.
for (col in score_cols) {
  base_col <- paste0(col, "_v1")
  if (base_col %in% names(wide)) {
    zero_baseline <- wide %>% filter(.data[[base_col]] == 0)
    if (nrow(zero_baseline) > 0) {
      warning(paste0(col, ": baseline score of 0 found for record(s) ",
                      paste(zero_baseline$record_id, collapse = ", "),
                      " - % change is undefined (division by zero) for these."))
    }
  }
}

pct_change_cols <- character(0)
for (col in score_cols) {
  base_col <- paste0(col, "_v1")
  fu_col   <- paste0(col, "_v6")
  pct_col  <- paste0(col, "_pctchange_v1_to_v6")
  if (all(c(base_col, fu_col) %in% names(wide))) {
    wide[[pct_col]] <- (wide[[fu_col]] - wide[[base_col]]) / wide[[base_col]] * 100
  } else {
    wide[[pct_col]] <- NA_real_
  }
  pct_change_cols <- c(pct_change_cols, pct_col)
}

# Composite: average of available % changes. Per PI direction: require ALL 3
# subtests evaluable (no missing data tolerated), given the reduced subtest
# count - stricter than the original protocol's proportional "4 of 5" rule.
wide <- wide %>%
  rowwise() %>%
  mutate(
    n_valid_subtests = sum(!is.na(c_across(all_of(pct_change_cols)))),
    composite_pct_change_v1_to_v6 = if (n_valid_subtests >= 3) {
      mean(c_across(all_of(pct_change_cols)), na.rm = TRUE)
    } else {
      NA_real_
    },
    unevaluable = n_valid_subtests < 3,
    adverse_cognitive_event = if (unevaluable) {
      NA
    } else {
      composite_pct_change_v1_to_v6 <= -10
    }
  ) %>%
  ungroup()

write.csv(wide, "adverse_event_pct_change_v1_to_v6.csv", row.names = FALSE)

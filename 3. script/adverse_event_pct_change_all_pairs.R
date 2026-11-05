library(dplyr)
library(tidyr)
library(stringr)

# --- 1. Read the export ------------------------------------------------------
df <- read.csv("~/Downloads/NIHTB_corrected_raw_3008.csv", stringsAsFactors = FALSE)

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

# --- 4. Pivot ALL visits wide (not just v1/v6) -------------------------------
# Needed once, up front, so each of the three visit-pair blocks below can pull
# from the same wide table rather than re-filtering/re-pivoting per pair.
wide <- clean %>%
  pivot_wider(id_cols = record_id, names_from = visit, values_from = all_of(score_cols))

# Flag zero/near-zero baselines across ALL visits that could serve as a
# "from" timepoint (v1 and v6, since v6 is the baseline for the v6->v12 pair).
# % change is undefined (division by zero) or explosively unstable otherwise.
for (from_visit in c("v1", "v6")) {
  for (col in score_cols) {
    base_col <- paste0(col, "_", from_visit)
    if (base_col %in% names(wide)) {
      zero_baseline <- wide %>% filter(.data[[base_col]] == 0)
      if (nrow(zero_baseline) > 0) {
        warning(paste0(col, " at ", from_visit, ": score of 0 found for record(s) ",
                        paste(zero_baseline$record_id, collapse = ", "),
                        " - % change from this visit is undefined for these."))
      }
    }
  }
}

# --- 5. Compute % change, composite, and AE flag for one visit pair ---------
# Shared logic for all three pairs (v1->v6, v6->v12, v1->v12), so the rule is
# defined once and applied identically rather than copy-pasted three times
# with room for the rules to silently drift apart between pairs.
add_pct_change_block <- function(wide_df, from_visit, to_visit, score_cols) {
  suffix <- paste0(from_visit, "_to_", to_visit)
  pct_change_cols <- character(0)

  for (col in score_cols) {
    base_col <- paste0(col, "_", from_visit)
    fu_col   <- paste0(col, "_", to_visit)
    pct_col  <- paste0(col, "_pctchange_", suffix)

    if (all(c(base_col, fu_col) %in% names(wide_df))) {
      wide_df[[pct_col]] <- (wide_df[[fu_col]] - wide_df[[base_col]]) / wide_df[[base_col]] * 100
    } else {
      wide_df[[pct_col]] <- NA_real_
    }
    pct_change_cols <- c(pct_change_cols, pct_col)
  }

  n_valid_col  <- paste0("n_valid_subtests_", suffix)
  comp_col     <- paste0("composite_pct_change_", suffix)
  uneval_col   <- paste0("unevaluable_", suffix)
  ae_col       <- paste0("adverse_cognitive_event_", suffix)

  wide_df <- wide_df %>%
    rowwise() %>%
    mutate(
      !!n_valid_col := sum(!is.na(c_across(all_of(pct_change_cols)))),
      !!comp_col := if (get(n_valid_col) >= 3) {
        mean(c_across(all_of(pct_change_cols)), na.rm = TRUE)
      } else {
        NA_real_
      },
      !!uneval_col := get(n_valid_col) < 3,
      !!ae_col := if (get(uneval_col)) {
        NA
      } else {
        get(comp_col) <= -10
      }
    ) %>%
    ungroup()

  wide_df
}

# --- 6. Apply to all three visit pairs ---------------------------------------
# Missing-data rule (per PI direction): require ALL 3 subtests evaluable at
# BOTH ends of a given pair - applied identically to each pair. A record can
# be evaluable for v1->v6 but unevaluable for v6->v12 (or vice versa) if data
# is missing at only one of the relevant visits; that's expected, not a bug.
wide <- add_pct_change_block(wide, "v1", "v6",  score_cols)
wide <- add_pct_change_block(wide, "v6", "v12", score_cols)
wide <- add_pct_change_block(wide, "v1", "v12", score_cols)

write.csv(wide, "adverse_event_pct_change_all_pairs_3008.csv", row.names = FALSE)

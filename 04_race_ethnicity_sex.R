## =====================================================================
## Race, ethnicity and sex for Table 1 / Table 2 (Annals policy)
## Stand-alone: uses 07_analytic_dataset.csv from the main analysis and
## visits.csv from MC-MED. No re-run of the main analysis is needed.
## Output: 20_race_ethnicity_sex.csv (trial-visit level, like Tables 1-2)
##         20b_race_ethnicity_sex_patients.csv (patient level, first eligible visit)
## =====================================================================
suppressPackageStartupMessages(library(data.table))

home <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else path.expand("~")
data_dir <- Sys.getenv("MCMED_DATA_DIR", unset = file.path(home, "Desktop", "mimiciv", "MC_MED"))
out_dir  <- Sys.getenv("MCMED_OUT_DIR", unset = file.path(home, "Desktop", "mimiciv", "analysis_dap_pp_bolus_seqtrials"))

d   <- fread(file.path(out_dir, "07_analytic_dataset.csv"), colClasses = list(character = "id"))
vc  <- names(fread(file.path(data_dir, "visits.csv"), nrows = 0))
pid <- intersect(c("MRN", "Patient_id", "PatientID"), vc)[1]
need <- c("CSN", pid, "Gender", "Race", "Ethnicity")
miss <- setdiff(need[!is.na(need)], vc)
if (length(miss)) stop("Columns not found in visits.csv: ", paste(miss, collapse = ", "))
vis <- fread(file.path(data_dir, "visits.csv"), select = need[!is.na(need)],
             colClasses = list(character = "CSN"))
setnames(vis, "CSN", "id")
if (!is.na(pid)) setnames(vis, pid, "pid") else vis[, pid := id]

# Categories exactly as recorded; missing, "Unknown" and "Declines to State" pooled
clean <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | grepl("unknown|declin|unavailable|^u$", x, ignore.case = TRUE)] <-
    "Unknown or declined"
  x
}
vis[, `:=`(sex = fifelse(Gender == "F", "Female", fifelse(Gender == "M", "Male", "Unknown")),
           race = clean(Race), ethnicity = clean(Ethnicity))]

message("Race categories in MC-MED:");      print(vis[id %in% d$id, .N, by = race][order(-N)])
message("Ethnicity categories in MC-MED:"); print(vis[id %in% d$id, .N, by = ethnicity][order(-N)])

x <- merge(d[, .(id, k, bolus, pheno)], vis[, .(id, pid, sex, race, ethnicity)], by = "id", all.x = TRUE)

tab <- function(dt, group) {
  rbind(
    dt[, .(n = .N), by = .(category = sex)][, variable := "Sex"],
    dt[, .(n = .N), by = .(category = race)][, variable := "Race"],
    dt[, .(n = .N), by = .(category = ethnicity)][, variable := "Ethnicity"]
  )[, `:=`(group = group, value = sprintf("%s (%.1f)", format(n, big.mark = ","), 100 * n / sum(n))), by = variable][]
}
groups <- list(Overall = x, Bolus = x[bolus == 1], `No bolus` = x[bolus == 0])
for (p in unique(na.omit(x$pheno))) groups[[p]] <- x[pheno == p]
long <- rbindlist(Map(tab, groups, names(groups)))
wide <- dcast(long, variable + category ~ group, value.var = "value", fill = "0 (0.0)")
setcolorder(wide, intersect(c("variable", "category", names(groups)), names(wide)))
wide[, variable := factor(variable, c("Sex", "Race", "Ethnicity"))]
setorder(wide, variable, category)
fwrite(wide, file.path(out_dir, "20_race_ethnicity_sex.csv"))

# patient level (each patient counted once, at the first eligible visit)
xp <- x[order(pid, id, k)][, .SD[1], by = pid]
fwrite(dcast(tab(xp, "Patients"), variable + category ~ group, value.var = "value"),
       file.path(out_dir, "20b_race_ethnicity_sex_patients.csv"))
print(wide)
message("Saved 20_race_ethnicity_sex.csv and 20b_race_ethnicity_sex_patients.csv in ", out_dir)

# =====================================================================
# MC-MED — Does the DAP / pulse-pressure phenotype at septic hypotension
# modify the effect of a crystalloid bolus?  Sequential target trials
# ---------------------------------------------------------------------
# Population : ED visits with infection (ICD-10 sepsis code or blood
#              culture + antibiotic within 24 h) and cuff hypotension
#              (MAP < 65 or SBP < 90) within 6 h of rooming
# Trials     : a new trial every 15 min from the first hypotensive reading
#              up to 3 h later. A visit is eligible for trial k if, at its
#              start: the latest cuff reading (previous 15 min) is
#              hypotensive; no norepinephrine started before trial start; no
#              crystalloid bolus in the previous 30 min; still in the ED.
#              Eligibility and confounders use ONLY information available at
#              trial start (revision 1: no conditioning on the grace window).
#              Visits can enter several trials (prior boluses allowed)
# Exposure   : crystalloid IV bolus started within 15 min of trial start
# Outcomes   : primary  - MAP at 60 min (median cuff MAP 45-75 min after start)
#              secondary - perfusion index at 60 min (ratio; trials with Perf at
#              start), norepinephrine within 6 h; MAP/Perf trajectories -30 to 120 min
# Extras     : number of prior boluses (0, 1, >=2) as pre-specified modifier;
#              diastolic shock index (HR/DAP) as continuous modifier; approximate
#              Bayesian analysis (normal priors); alluvial plot of phenotype transitions;
#              extra figures A (MAP change vs MAP at start) and B (probability of benefit)
# Modifiers  : phenotype on the trial-start reading (low pulse pressure =
#              DAP > 50 & PP < 40; vasodilatory = DAP <= 50 & PP >= 40;
#              mixed = both low) and continuous DAP and PP
# Analysis   : pooled over trials; stabilized IPW (truncated 1st-99th pct)
#              + weighted outcome model with exposure x modifier interaction,
#              g-computation within phenotype; cluster bootstrap resampling
#              patients (MRN) for 95% CIs
# Requires   : data.table, arrow, dplyr, lubridate, rlang, ggplot2, splines, ggalluvial
# =====================================================================


suppressPackageStartupMessages({
  library(data.table); library(arrow); library(dplyr); library(lubridate)
  library(rlang); library(splines); library(ggplot2)
})
set.seed(20261008)

## ---- 0. Configuration -------------------------------------------------
home <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else path.expand("~")
cfg <- list(
  data_dir = Sys.getenv("MCMED_DATA_DIR", unset = file.path(home, "Desktop", "mimiciv", "MC_MED")),
  out_dir  = Sys.getenv("MCMED_OUT_DIR", unset = file.path(home, "Desktop", "mimiciv", "analysis_dap_pp_bolus_seqtrials")),
  f_visits = "visits.csv", f_orders = "orders.csv", f_labs = "labs.csv", f_numerics = "numerics.csv",

  m    = list(sbp = "SBP", dbp = "DBP", map = "MAP", hr = "HR", rr = "RR", spo2 = "SpO2", perf = "Perf"),
  lims = list(sbp = c(40, 260), dbp = c(15, 160), map = c(25, 200), hr = c(20, 250),
              rr = c(4, 60), spo2 = c(50, 100), perf = c(0.01, 25)),

  horizon_min = 360,  si_window_h = 24,
  map_hypo = 65, sbp_hypo = 90, dap_low = 50, pp_low = 40, pp_min = 10,
  trial_step  = 15,    # a new trial starts every 15 min after index
  max_trial   = 180,   # last trial starts 3 h after index
  grace       = 15,    # exposure window: bolus started in [trial start, +15 min)
  washout     = 30,    # no crystalloid bolus in the 30 min before trial start
  out_win     = c(45, 75), # MAP at 60 min = median cuff MAP 45-75 min after trial start
  ne_horizon  = 360,   # norepinephrine within 6 h of trial start
  cov_win     = 10,    # minutes before trial start for monitor covariates
  n_boot = 300, n_boot_sens = 200, w_trunc = c(0.01, 0.99),
  mcid_map = 2,        # mm Hg: smallest clinically relevant MAP difference (also ROPE half-width)
  mcid_ne  = 0.05,
  mcid_perf = 0.10,    # relative change in perfusion index (10%)
  n_boot_traj = 200,   # bootstrap replicates for the trajectory figure     # risk difference: smallest relevant change in norepinephrine use

  rx_sepsis_icd = "A40|A41|R65\\.?2",
  rx_culture    = "blood culture|culture,? blood",
  rx_abx        = paste("cefepime|ceftriaxone|piperacillin|meropenem|vancomycin",
                        "ertapenem|imipenem|aztreonam|ceftazidime|cefazolin|cefuroxime",
                        "ampicillin|azithromycin|levofloxacin|ciprofloxacin",
                        "metronidazole|linezolid|daptomycin|gentamicin|clindamycin",
                        "doxycycline", sep = "|"),
  rx_ne         = "norepinephrine",
  rx_fluid      = "^(LR|NS|NORMOSOL-R|ELECTROLYTE-R|LACTATED RINGERS) IV BOLUS",
  rx_fluid_excl = "PREHYDRATION",
  rx_lactate    = "LACTATE|LACTIC",   # any lactate component (POC, venous, arterial...)
  rx_lactate_excl = "DEHYDROGENASE|LDH|CSF|PLEURAL|PERITONEAL|FLUID"
)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(file.path(cfg$data_dir, cfg$f_visits))) {
  hit <- list.files(cfg$data_dir, pattern = "^visits\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(hit)) stop("visits.csv not found under ", cfg$data_dir, call. = FALSE)
  cfg$data_dir <- dirname(hit[1])
}
fp   <- function(f) file.path(cfg$data_dir, f)
step <- function(txt) message(format(Sys.time(), "%H:%M:%S"), " — ", txt)
parse_dt <- function(x) {
  if (inherits(x, "POSIXct")) return(force_tz(x, "UTC"))
  suppressWarnings(parse_date_time(as.character(x),
    orders = c("Ymd HMS", "Ymd HM", "mdY HMS", "mdY HM"), tz = "UTC"))
}
mins <- function(a, b) as.numeric(difftime(a, b, units = "mins"))

## ---- 1. Visits and infection ------------------------------------------
step("1. Visits and infection")
vis_cols <- names(fread(fp(cfg$f_visits), nrows = 0))
pid_col  <- intersect(c("MRN", "Patient_id", "PatientID"), vis_cols)[1]   # patient identifier, if available
vis <- fread(fp(cfg$f_visits), select = c("CSN", if (!is.na(pid_col)) pid_col, "Age", "Gender", "Triage_acuity", "Means_of_arrival",
                                         "ED_dispo", "Dx_ICD10", "Arrival_time", "Roomed_time", "Departure_time"))
setnames(vis, "CSN", "id"); vis[, id := as.character(id)]
if (!is.na(pid_col)) { setnames(vis, pid_col, "pid"); vis[, pid := as.character(pid)] } else vis[, pid := id]
message("Bootstrap clusters: ", if (is.na(pid_col)) "visits (no patient identifier found)" else paste("patients (", pid_col, ")"))
vis[, t0 := parse_dt(Roomed_time)]; vis[is.na(t0), t0 := parse_dt(Arrival_time)]
vis[, dep_rel := mins(parse_dt(Departure_time), t0)]
vis[, sepsis_dx := grepl(cfg$rx_sepsis_icd, Dx_ICD10, ignore.case = TRUE)]

ord <- fread(fp(cfg$f_orders), select = c("CSN", "Order_time", "Procedure_name", "First_admin_time"))
setnames(ord, c("id", "o_time", "name", "a_time")); ord[, id := as.character(id)]
ord[, `:=`(o_time = parse_dt(o_time), a_time = parse_dt(a_time))]
cult <- ord[grepl(cfg$rx_culture, name, ignore.case = TRUE), .(t_cult = min(o_time, na.rm = TRUE)), by = id]
abx  <- ord[grepl(cfg$rx_abx,     name, ignore.case = TRUE), .(t_abx  = min(o_time, na.rm = TRUE)), by = id]
si   <- merge(cult, abx, by = "id")[abs(mins(t_abx, t_cult)) <= cfg$si_window_h * 60]
vis[, infection := sepsis_dx | id %in% si$id]
inf_ids <- vis[infection == TRUE & !is.na(t0), id]

## ---- 2. Administrations -------------------------------------------------
step("2. Fluids and norepinephrine")
adm <- ord[id %in% inf_ids & !is.na(a_time)]
rm(ord); gc(verbose = FALSE)
adm <- merge(adm, vis[, .(id, t0)], by = "id")
adm[, t_rel := mins(a_time, t0)]
ne  <- adm[grepl(cfg$rx_ne, name, ignore.case = TRUE), .(t_ne = min(t_rel, na.rm = TRUE)), by = id]
bol <- adm[grepl(cfg$rx_fluid, name, ignore.case = TRUE) & !grepl(cfg$rx_fluid_excl, name, ignore.case = TRUE),
           .(id, t_b = t_rel, name)]
bol[, bolus_type := fifelse(grepl("ML/KG", name, ignore.case = TRUE), "30 mL/kg",
                     fifelse(grepl("-\\s*\\d+\\s*ML", name, ignore.case = TRUE),
                             paste(sub(".*-\\s*(\\d+)\\s*ML.*", "\\1", toupper(name)), "mL"), "unspecified"))]

## ---- 3. Monitor and cuff numerics -------------------------------------
step("3. Numerics (this step takes several minutes)")
map_meas  <- rbindlist(lapply(names(cfg$m), function(k) data.table(measure = cfg$m[[k]], var = k)))
ids_keep  <- as.integer(inf_ids); meas_keep <- unique(map_meas$measure)
nm <- open_dataset(fp(cfg$f_numerics), format = "csv") |>
  filter(Measure %in% !!meas_keep) |>
  filter(CSN %in% !!ids_keep) |>
  select(CSN, Measure, Value, Time) |>
  collect() |> as.data.table()
setnames(nm, c("id", "measure", "value", "time"))
nm[, `:=`(id = as.character(id), value = as.numeric(value), time = parse_dt(time))]
nm <- merge(nm, map_meas, by = "measure")
nm <- merge(nm, vis[, .(id, t0)], by = "id")
nm[, t_rel := mins(time, t0)]
nm <- nm[!is.na(value) & t_rel >= -60 & t_rel <= cfg$horizon_min + 120]
for (v in names(cfg$lims)) nm <- nm[!(var == v & (value < cfg$lims[[v]][1] | value > cfg$lims[[v]][2]))]

## ---- 4. Cuff readings and index hypotension ---------------------------
step("4. Cuff readings and index hypotension")
nbp <- dcast(nm[var %in% c("sbp", "dbp", "map"), .(id, tm = round(t_rel), var, value)],
             id + tm ~ var, value.var = "value", fun.aggregate = median)
nbp <- nbp[!is.na(sbp) & !is.na(dbp) & sbp - dbp >= cfg$pp_min]
nbp[, map := fifelse(!is.na(map), map, (sbp + 2 * dbp) / 3)]

idx <- nbp[tm >= 0 & tm <= cfg$horizon_min & (map < cfg$map_hypo | sbp < cfg$sbp_hypo)][order(tm), .SD[1], by = id]
idx <- idx[, .(id, t_idx = tm)]

## ---- 5. Sequential trials -------------------------------------------------
step("5. Sequential trials")
base <- merge(idx[, .(id, t_idx)],
              vis[, .(id, pid, Age, Gender, Triage_acuity, Means_of_arrival, ED_dispo, dep_rel)], by = "id")
base <- merge(base, ne, by = "id", all.x = TRUE)
tr <- base[, .(k = seq(0, cfg$max_trial, by = cfg$trial_step)), by = id]
tr <- merge(tr, base, by = "id")
tr[, t_start := t_idx + k]
tr <- tr[t_start <= cfg$horizon_min]
flow <- list(c("Infection visits", length(inf_ids)),
             c("Visits with cuff hypotension within 6 h", nrow(base)),
             c("Candidate trial-visits (every 15 min, up to 3 h after index)", nrow(tr)))

# revision 1: eligibility uses only information available at trial start
tr <- tr[is.na(t_ne)    | t_ne    > t_start]; flow[[4]] <- c("  no norepinephrine started before trial start", nrow(tr))
tr <- tr[is.na(dep_rel) | dep_rel > t_start]; flow[[5]] <- c("  still in ED at trial start", nrow(tr))

# latest cuff reading in the 15 min before trial start, must be hypotensive
cur <- nbp[tr[, .(id, k, t_start)], on = "id", allow.cartesian = TRUE, nomatch = 0][
  tm <= t_start & tm > t_start - cfg$trial_step]
cur <- cur[order(-tm), .SD[1], by = .(id, k)][, .(id, k, t_cur = tm, sbp_cur = sbp, dap = dbp, map_cur = map, pp = sbp - dbp)]
tr  <- merge(tr, cur, by = c("id", "k"));                 flow[[6]] <- c("  cuff reading in previous 15 min", nrow(tr))
tr  <- tr[map_cur < cfg$map_hypo | sbp_cur < cfg$sbp_hypo]; flow[[7]] <- c("  currently hypotensive", nrow(tr))

# bolus history and exposure
bh <- bol[tr[, .(id, k, t_start)], on = "id", allow.cartesian = TRUE, nomatch = 0]
recent <- bh[t_b < t_start & t_b >= t_start - cfg$washout, .(id, k)][, recent := TRUE]
tr <- merge(tr, unique(recent), by = c("id", "k"), all.x = TRUE)[is.na(recent)][, recent := NULL]
flow[[8]] <- c("  no bolus in previous 30 min (eligible trial-visits)", nrow(tr))
nprior <- bh[t_b < t_start, .(n_prior = .N), by = .(id, k)]
expo   <- bh[t_b >= t_start & t_b < t_start + cfg$grace][order(t_b), .SD[1], by = .(id, k)][, .(id, k, bolus_type)]
late   <- bh[t_b >= t_start + cfg$grace & t_b < t_start + 60, .(id, k)][, late_b := TRUE]
tr <- merge(tr, nprior, by = c("id", "k"), all.x = TRUE); tr[is.na(n_prior), n_prior := 0L]
tr <- merge(tr, expo,   by = c("id", "k"), all.x = TRUE); tr[, bolus := as.integer(!is.na(bolus_type))]
tr <- merge(tr, unique(late), by = c("id", "k"), all.x = TRUE)
tr[, ctrl_late_bolus := bolus == 0 & !is.na(late_b)]; tr[, late_b := NULL]

# phenotype on the trial-start reading
tr[, pheno := fcase(dap <= cfg$dap_low & pp >= cfg$pp_low, "Vasodilatory",
                    dap >  cfg$dap_low & pp <  cfg$pp_low, "Low pulse pressure",
                    dap <= cfg$dap_low & pp <  cfg$pp_low, "Mixed",
                    default = "Other")]
flow[[9]] <- c("  phenotype other (DAP > 50 and PP >= 40), excluded", tr[pheno == "Other", .N])
tr <- tr[pheno != "Other"]
tr[, pheno := factor(pheno, levels = c("Low pulse pressure", "Vasodilatory", "Mixed"))]

## ---- 6. Covariates and outcomes -------------------------------------------
step("6. Covariates and outcomes")
near_tr <- function(v) {
  s <- nm[var == v, .(id, t_rel, value)]
  q <- tr[, .(id, lo = t_start - cfg$cov_win, hi = t_start)]
  s[q, on = .(id, t_rel >= lo, t_rel <= hi), median(value), by = .EACHI]$V1
}
for (v in c("hr", "rr", "spo2", "perf")) tr[, (v) := near_tr(v)]

pv <- nbp[tr[, .(id, k, t_cur)], on = "id", allow.cartesian = TRUE, nomatch = 0][tm < t_cur & tm >= t_cur - 60]
pv <- pv[order(-tm), .SD[1], by = .(id, k)][, .(id, k, map_prev = map)]
tr <- merge(tr, pv, by = c("id", "k"), all.x = TRUE)

# lactate: latest specimen ORDERED up to trial start (value used even if resulted later; revision 1)
lab <- fread(fp(cfg$f_labs), select = c("CSN", "Order_time", "Component_name", "Component_value"))
setnames(lab, c("id", "o_time", "comp", "val")); lab[, id := as.character(id)]
lab <- lab[id %in% unique(tr$id) & grepl(cfg$rx_lactate, comp, ignore.case = TRUE) &
           !grepl(cfg$rx_lactate_excl, comp, ignore.case = TRUE)]
fwrite(lab[, .N, by = comp][order(-N)], file.path(cfg$out_dir, "00_lactate_components.csv"))   # check
lab[, `:=`(val = suppressWarnings(as.numeric(gsub("[^0-9.]", "", val))), o_time = parse_dt(o_time))]
lab <- merge(lab, vis[, .(id, t0)], by = "id")[, t_rel := mins(o_time, t0)][!is.na(val)]
lt  <- lab[tr[, .(id, k, t_start)], on = "id", allow.cartesian = TRUE, nomatch = 0][t_rel <= t_start]
lt  <- lt[order(-t_rel), .SD[1], by = .(id, k)][, .(id, k, lactate = val)]
tr  <- merge(tr, lt, by = c("id", "k"), all.x = TRUE)
rm(lab); gc(verbose = FALSE)

om <- nbp[tr[, .(id, k, t_start)], on = "id", allow.cartesian = TRUE, nomatch = 0][
  tm - t_start >= cfg$out_win[1] & tm - t_start <= cfg$out_win[2], .(map60 = median(map)), by = .(id, k)]
tr <- merge(tr, om, by = c("id", "k"), all.x = TRUE)
op <- nm[var == "perf", .(id, t_rel, value)][tr[, .(id, k, t_start)], on = "id", allow.cartesian = TRUE, nomatch = 0][
  t_rel - t_start >= cfg$out_win[1] & t_rel - t_start <= cfg$out_win[2], .(perf60 = median(value)), by = .(id, k)]
tr <- merge(tr, op, by = c("id", "k"), all.x = TRUE)
tr[, lperf60 := log(perf60)]
tr[, ne6h := as.integer(!is.na(t_ne) & t_ne <= t_start + cfg$ne_horizon)]
tr[, icu  := as.integer(ED_dispo %in% "ICU")]

tr[, `:=`(age = Age, female = as.integer(Gender == "F"),
          acuity_hi = as.integer(suppressWarnings(as.integer(substr(Triage_acuity, 1, 1))) <= 2),
          ems = as.integer(grepl("EMS|ambulance", Means_of_arrival, ignore.case = TRUE)),
          prior_bolus = as.integer(n_prior > 0), n_prior_c = pmin(n_prior, 3L))]
tr[is.na(acuity_hi), acuity_hi := 0L]
tr[, d_prev := map_prev - map_cur]
tr[, `:=`(hr_raw = hr, perf_raw = perf, lactate_raw = lactate)]
imp <- function(x) { m <- is.na(x); x[m] <- median(x, na.rm = TRUE); list(x, as.integer(m)) }
for (v in c("hr", "rr", "spo2", "perf", "lactate", "d_prev", "age")) tr[, c(v, paste0(v, "_miss")) := imp(get(v))]
tr[, lperf := log(perf)]
tr[, dsi := hr / dap]
tr[, nprior3 := factor(fcase(n_prior == 0, "0", n_prior == 1, "1", default = ">=2"), levels = c("0", "1", ">=2"))]          # diastolic shock index (HR / DAP), HR median-imputed if missing
# time-varying information DURING the grace window (used only by the clone-censor-weight sensitivity analysis)
gs <- seq(0, cfg$grace - 5, by = 5)                                  # starts of 5-min intervals: 0, 5, 10
fb <- bh[t_b >= t_start & t_b < t_start + cfg$grace, .(b_rel = min(t_b - t_start)), by = .(id, k)]
tr <- merge(tr, fb, by = c("id", "k"), all.x = TRUE)
for (g in gs) {
  tg <- tr[, .(id, k, tg = t_start + g, map_cur)]
  rg <- nbp[tg, on = "id", allow.cartesian = TRUE, nomatch = 0][tm <= tg & tm > tg - cfg$trial_step]
  rg <- rg[order(-tm), .SD[1], by = .(id, k)][, .(id, k, mg = map)]
  tr <- merge(tr, rg, by = c("id", "k"), all.x = TRUE)
  tr[is.na(mg), mg := map_cur]
  setnames(tr, "mg", paste0("map_g", g))
  tr[, paste0("ne_g", g) := as.integer(!is.na(t_ne) & t_ne <= t_start + g)]
  tr[, paste0("in_ed_g", g) := as.integer(is.na(dep_rel) | dep_rel > t_start + g)]
}
d <- tr
flow[[10]] <- c("Analytic trial-visits", nrow(d))
flow[[11]] <- c("  distinct visits", uniqueN(d$id))
flow[[17]] <- c("  distinct patients", uniqueN(d$pid))
flow[[12]] <- c("  trial-visits with bolus (exposed)", d[bolus == 1, .N])
flow[[13]] <- c("  distinct visits ever exposed", d[bolus == 1, uniqueN(id)])
flow[[14]] <- c("  trial-visits with MAP at 60 min", d[!is.na(map60), .N])
flow[[15]] <- c("  trial-visits with prior bolus", d[prior_bolus == 1, .N])
flow[[16]] <- c("  trial-visits with perfusion index at start and at 60 min", d[perf_miss == 0 & !is.na(lperf60), .N])

## ---- 7. Models ---------------------------------------------------------------
step("7. Models and cluster bootstrap (several minutes)")
K <- list(dap = quantile(d$dap, c(.33, .67)), pp = quantile(d$pp, c(.33, .67)), dsi = quantile(d$dsi, c(.33, .67)))
B <- list(dap = range(d$dap), pp = range(d$pp), dsi = range(d$dsi))
sp   <- "ns(dap, knots = K$dap, Boundary.knots = B$dap) + ns(pp, knots = K$pp, Boundary.knots = B$pp)"
covs <- paste("age + female + acuity_hi + ems + ns(t_start, df = 3) + k + n_prior_c",
              "+ hr + hr_miss + rr + rr_miss + spo2 + spo2_miss",
              "+ lperf + perf_miss + lactate + lactate_miss + d_prev + d_prev_miss")
f_ps    <- as.formula(paste("bolus ~", sp, "+", covs))
f_map   <- as.formula(paste("map60 ~ bolus * pheno +", sp, "+", covs))
f_ne    <- as.formula(paste("ne6h ~ bolus * pheno +", sp, "+", covs))
f_map_c <- as.formula(paste("map60 ~ bolus * (", sp, ") +", covs))
f_ne_c  <- as.formula(paste("ne6h ~ bolus * (", sp, ") +", covs))
f_perf  <- as.formula(paste("lperf60 ~ bolus * pheno +", sp, "+", covs))
# pre-specified modifier: number of crystalloid boluses already received (0, 1, >=2)
covs_np  <- sub("\\+ n_prior_c", "", covs)
f_map_n  <- as.formula(paste("map60 ~ bolus * nprior3 + pheno +", sp, "+", covs_np))
f_ne_n   <- as.formula(paste("ne6h ~ bolus * nprior3 + pheno +", sp, "+", covs_np))
f_perf_n <- as.formula(paste("lperf60 ~ bolus * nprior3 + pheno +", sp, "+", covs_np))
sdsi    <- "ns(dsi, knots = K$dsi, Boundary.knots = B$dsi)"
f_map_d <- as.formula(paste("map60 ~ bolus *", sdsi, "+", sp, "+", covs))
f_ne_d  <- as.formula(paste("ne6h ~ bolus *",  sdsi, "+", sp, "+", covs))
grid <- list(dap = seq(quantile(d$dap, .05), quantile(d$dap, .95), length.out = 20),
             pp  = seq(quantile(d$pp,  .05), quantile(d$pp,  .95), length.out = 20),
             dsi = seq(quantile(d$dsi, .05), quantile(d$dsi, .95), length.out = 20))

ipw <- function(x) {
  ps <- fitted(glm(f_ps, family = binomial(), data = x))
  p1 <- mean(x$bolus)
  w  <- ifelse(x$bolus == 1, p1 / ps, (1 - p1) / (1 - ps))
  q  <- quantile(w, cfg$w_trunc); pmin(pmax(w, q[1]), q[2])
}
# Clone-censor-weight for the 15-min grace window. Clones of the "bolus within 15 min" arm are
# uncensored only if a bolus starts in the window (weight 1 / P(bolus in window | baseline));
# clones of the "no bolus" arm are censored when a bolus starts, with weights from a pooled
# logistic model of the bolus hazard in 5-min intervals given baseline AND time-varying
# information (norepinephrine started, latest cuff MAP). Uncensored clones are disjoint across arms.
f_haz <- as.formula(paste("ev ~", sp, "+", covs, "+ factor(g) + ne_on + map_g"))
ipw_tv <- function(x) {
  x <- copy(x)[, row := .I]
  lg <- rbindlist(lapply(gs, function(g) x[, .(row, g = g, ne_on = get(paste0("ne_g", g)), map_g = get(paste0("map_g", g)),
                                               at_risk = get(paste0("in_ed_g", g)) == 1 & (is.na(b_rel) | b_rel >= g),
                                               ev = as.integer(!is.na(b_rel) & b_rel >= g & b_rel < g + 5))]))
  lg <- lg[at_risk == TRUE]
  cols <- setdiff(intersect(names(x), all.vars(f_ps)), "bolus")
  lg <- merge(lg, x[, c("row", cols), with = FALSE], by = "row")
  lg[, h := fitted(glm(f_haz, family = binomial(), data = lg))]
  # "no bolus" clones: censored at bolus start -> P(uncensored) = product of (1 - hazard) over the intervals at risk
  p_nob  <- lg[, .(p0 = prod(1 - h)), by = row]
  # "bolus" clones: censored at the END of the window if still untreated. Only clones untreated at the start of the
  # last interval are at risk of this censoring; their P(uncensored) = hazard in the last interval. Clones treated in
  # earlier intervals cannot be censored (weight 1). Clones leaving the ED untreated are censored at departure (unweighted).
  last_g <- max(gs)
  h_last <- lg[g == last_g, .(row, h_last = h)]
  x <- merge(x, p_nob, by = "row", all.x = TRUE); x <- merge(x, h_last, by = "row", all.x = TRUE)
  x[is.na(p0), p0 := 1]; x[is.na(h_last), h_last := 1]; setorder(x, row)   # h_last missing: not at risk (e.g. left ED)
  w <- ifelse(x$bolus == 1, ifelse(x$b_rel >= last_g, 1 / pmax(x$h_last, 0.01), 1), 1 / pmax(x$p0, 0.01))
  w[x$bolus == 1] <- w[x$bolus == 1] / mean(w[x$bolus == 1])     # normalized within arm (mean 1)
  w[x$bolus == 0] <- w[x$bolus == 0] / mean(w[x$bolus == 0])
  q  <- quantile(w, cfg$w_trunc); pmin(pmax(w, q[1]), q[2])
}
# Phenotype-stratified propensity model (post hoc): separate models within each phenotype, splines with 3 df
f_ps_ph <- as.formula(paste("bolus ~ ns(dap, df = 3) + ns(pp, df = 3) +", covs))
ipw_ph <- function(x) {
  w <- numeric(nrow(x))
  for (p in levels(x$pheno)) {
    i <- which(x$pheno == p); xi <- x[i]
    ps <- fitted(glm(f_ps_ph, family = binomial(), data = xi)); p1 <- mean(xi$bolus)
    w[i] <- ifelse(xi$bolus == 1, p1 / ps, (1 - p1) / (1 - ps))
  }
  q <- quantile(w, cfg$w_trunc); pmin(pmax(w, q[1]), q[2])
}
gdiff <- function(model, x) {
  x1 <- copy(x)[, bolus := 1L]; x0 <- copy(x)[, bolus := 0L]
  predict(model, x1, type = "response") - predict(model, x0, type = "response")
}
f_obs_map  <- as.formula(paste("obs ~ bolus * pheno +", sp, "+", covs))
est <- function(x, weighted = TRUE, curves = TRUE, ipcw = FALSE, wfun = ipw) {
  x <- copy(x); x[, w := if (weighted) wfun(x) else 1]
  if (ipcw) {                  # inverse probability of having MAP at 60 min (post hoc)
    x[, obs := as.integer(!is.na(map60))]
    po <- fitted(glm(f_obs_map, family = binomial(), data = x))
    x[, w := w / pmax(po, 0.05)]
  }
  xm <- x[!is.na(map60)]
  m1 <- lm(f_map, data = xm, weights = w)
  m2 <- glm(f_ne, family = quasibinomial(), data = x, weights = w)
  xp <- x[perf_miss == 0 & !is.na(lperf60)]
  m5 <- lm(f_perf, data = xp, weights = w)
  em <- gdiff(m1, xm); en <- gdiff(m2, x); ep <- gdiff(m5, xp)   # ep: difference in log Perf (log ratio)
  r <- c(map_Overall = mean(em), ne_Overall = mean(en), perf_Overall = mean(ep))
  for (p in levels(x$pheno)) {
    r[paste0("map_", p)] <- mean(em[xm$pheno == p]); r[paste0("ne_", p)] <- mean(en[x$pheno == p])
    r[paste0("perf_", p)] <- mean(ep[xp$pheno == p])
  }
  for (o in c("map", "ne", "perf")) for (p in c("Vasodilatory", "Mixed"))
    r[paste0(o, "_contrast_", p, " vs Low pulse pressure")] <- r[paste0(o, "_", p)] - r[paste0(o, "_Low pulse pressure")]
  if (uniqueN(x$nprior3) == 3) {        # skipped in subsets where the modifier does not vary
    nm_lab <- c("0" = "Prior boluses: 0", "1" = "Prior boluses: 1", ">=2" = "Prior boluses: >=2")
    mods <- list(map  = list(lm(f_map_n, data = xm, weights = w), xm),
                 ne   = list(glm(f_ne_n, family = quasibinomial(), data = x, weights = w), x),
                 perf = list(lm(f_perf_n, data = xp, weights = w), xp))
    for (o in names(mods)) {
      gd <- gdiff(mods[[o]][[1]], mods[[o]][[2]]); lv <- mods[[o]][[2]]$nprior3
      for (l in names(nm_lab)) r[paste0(o, "_", nm_lab[l])] <- mean(gd[lv == l])
      for (l in c("1", ">=2")) r[paste0(o, "_contrast_", nm_lab[l], " vs 0")] <- mean(gd[lv == l]) - mean(gd[lv == "0"])
    }
  }
  if (curves) {
    m3 <- lm(f_map_c, data = xm, weights = w)
    m4 <- glm(f_ne_c, family = quasibinomial(), data = x, weights = w)
    m6 <- lm(f_map_d, data = xm, weights = w)
    m7 <- glm(f_ne_d, family = quasibinomial(), data = x, weights = w)
    for (g in seq_along(grid$dsi)) {
      xm_g <- copy(xm)[, dsi := grid$dsi[g]]; x_g <- copy(x)[, dsi := grid$dsi[g]]
      r[paste0("curve_map_dsi_", g)] <- mean(gdiff(m6, xm_g))
      r[paste0("curve_ne_dsi_",  g)] <- mean(gdiff(m7, x_g))
    }
    for (v in c("dap", "pp")) for (g in seq_along(grid[[v]])) {
      xm_g <- copy(xm)[, (v) := grid[[v]][g]]; x_g <- copy(x)[, (v) := grid[[v]][g]]
      r[paste0("curve_map_", v, "_", g)] <- mean(gdiff(m3, xm_g))
      r[paste0("curve_ne_",  v, "_", g)] <- mean(gdiff(m4, x_g))
    }
  }
  r
}
cluster_sample <- function(x) {
  s <- data.table(pid = sample(unique(x$pid), replace = TRUE))[, boot_id := .I]   # resample PATIENTS
  x[s, on = "pid", allow.cartesian = TRUE]
}
boot_ci <- function(x, nb, ...) {
  pt <- est(x, ...)
  bs <- t(vapply(seq_len(nb), function(b) {
    out <- tryCatch(est(cluster_sample(x), ...), error = function(e) NULL)
    if (is.null(out)) rep(NA_real_, length(pt)) else out[names(pt)]
  }, numeric(length(pt))))
  colnames(bs) <- names(pt)
  boot_draws[[length(boot_draws) + 1]] <<- bs          # kept for the Bayesian analysis
  data.table(term = names(pt), est = pt,
             se = apply(bs, 2, sd, na.rm = TRUE),
             lo = apply(bs, 2, quantile, .025, na.rm = TRUE),
             hi = apply(bs, 2, quantile, .975, na.rm = TRUE))
}
boot_draws <- list()

prim <- boot_ci(d, cfg$n_boot)
fwrite(as.data.table(boot_draws[[1]][, grep("^curve_", colnames(boot_draws[[1]]), invert = TRUE), drop = FALSE]),
       file.path(cfg$out_dir, "08_bootstrap_draws_primary.csv"))
step("   sensitivity analyses")
first_trial <- d[order(k), .SD[1], by = id]
s1 <- boot_ci(d[ctrl_late_bolus == FALSE], cfg$n_boot_sens, curves = FALSE)  # controls without bolus up to 60 min
s2 <- boot_ci(first_trial,                 cfg$n_boot_sens, curves = FALSE)  # first eligible trial per visit
s3 <- boot_ci(d[map_cur < cfg$map_hypo],   cfg$n_boot_sens, curves = FALSE)  # MAP < 65 at trial start
s4 <- boot_ci(d[prior_bolus == 0],         cfg$n_boot_sens, curves = FALSE)  # no prior bolus
s5 <- boot_ci(d, cfg$n_boot_sens, weighted = FALSE, curves = FALSE)          # outcome model without weights
s6 <- boot_ci(d, cfg$n_boot_sens, curves = FALSE, ipcw = TRUE)               # post hoc: IPCW for missing MAP at 60 min
# post hoc: dose (ordered volume). Exposed = boluses of that size only; trials with other boluses dropped
large_vol <- c("1000 mL", "2000 mL", "30 mL/kg"); small_vol <- c("250 mL", "500 mL")
s7 <- boot_ci(d[bolus == 0 | bolus_type %in% large_vol], cfg$n_boot_sens, curves = FALSE)
s8 <- boot_ci(d[bolus == 0 | bolus_type %in% small_vol], cfg$n_boot_sens, curves = FALSE)
s9 <- boot_ci(d, cfg$n_boot_sens, curves = FALSE, wfun = ipw_tv)            # revision: clone-censor-weight for grace window
# norepinephrine within 6 h is not estimable this way: events can occur before artificial censoring at the end of the window
s9[grepl("^ne_", term), c("est", "se", "lo", "hi") := NA_real_]
s10 <- boot_ci(d, cfg$n_boot_sens, curves = FALSE, wfun = ipw_ph)          # post hoc: phenotype-stratified propensity model

## ---- 8. Balance and descriptive tables ----------------------------------------
step("8. Tables")
d[, w := ipw(d)]
bal_vars <- c("age", "female", "acuity_hi", "ems", "t_start", "k", "n_prior_c", "dap", "pp", "map_cur", "hr", "rr",
              "spo2", "lperf", "perf_miss", "lactate", "lactate_miss", "d_prev")
smd <- function(x, t, w) {
  m1 <- weighted.mean(x[t == 1], w[t == 1]); m0 <- weighted.mean(x[t == 0], w[t == 0])
  v1 <- sum(w[t == 1] * (x[t == 1] - m1)^2) / sum(w[t == 1])
  v0 <- sum(w[t == 0] * (x[t == 0] - m0)^2) / sum(w[t == 0])
  (m1 - m0) / sqrt((v1 + v0) / 2)
}
balance <- rbindlist(lapply(bal_vars, function(v) data.table(
  variable = v, smd_unweighted = smd(d[[v]], d$bolus, rep(1, nrow(d))), smd_weighted = smd(d[[v]], d$bolus, d$w))))

mq <- function(x) sprintf("%.1f (%.1f-%.1f)", median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
np <- function(x) sprintf("%d (%.1f)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
desc <- function(x) data.table(
  characteristic = c("Trial-visits", "Distinct visits", "Age, years", "Female sex", "ESI 1-2", "EMS arrival",
                     "Minutes rooming to trial start", "Minutes index to trial start", "Prior crystalloid bolus",
                     "SBP at trial start, mm Hg", "DAP at trial start, mm Hg", "Pulse pressure at trial start, mm Hg",
                     "MAP at trial start, mm Hg", "Heart rate, bpm", "Perfusion index", "Perfusion index missing",
                     "Lactate, mmol/L", "Lactate missing",
                     "Phenotype: low pulse pressure", "Phenotype: vasodilatory", "Phenotype: mixed",
                     "MAP at 60 min, mm Hg", "Norepinephrine within 6 h", "ICU admission"),
  value = c(nrow(x), uniqueN(x$id), mq(x$Age), np(x$female), np(x$acuity_hi), np(x$ems),
            mq(x$t_start), mq(x$k), np(x$prior_bolus),
            mq(x$sbp_cur), mq(x$dap), mq(x$pp), mq(x$map_cur), mq(x$hr_raw), mq(x$perf_raw), np(x$perf_miss),
            mq(x$lactate_raw), np(x$lactate_miss),
            np(x$pheno == "Low pulse pressure"), np(x$pheno == "Vasodilatory"), np(x$pheno == "Mixed"),
            mq(x$map60), np(x$ne6h), np(x$icu)))
t1 <- Reduce(function(a, b) merge(a, b, by = "characteristic", sort = FALSE),
             list(setnames(desc(d), "value", "Overall"), setnames(desc(d[bolus == 1]), "value", "Bolus"),
                  setnames(desc(d[bolus == 0]), "value", "No bolus")))
t1_ph <- Reduce(function(a, b) merge(a, b, by = "characteristic", sort = FALSE),
                lapply(levels(d$pheno), function(p) setnames(desc(d[pheno == p]), "value", p)))

## ---- 9. Outputs ----------------------------------------------------------------
step("9. Outputs and figures")
eff <- rbind(prim[!grepl("^curve_", term)][, analysis := "Primary (sequential trials, IPW + outcome model)"],
             s1[, analysis := "Controls without bolus up to 60 min"],
             s2[, analysis := "First eligible trial per visit"],
             s3[, analysis := "MAP < 65 mm Hg at trial start"],
             s4[, analysis := "No prior bolus"],
             s5[, analysis := "Outcome model only (unweighted)"],
             s6[, analysis := "Post hoc: inverse probability of MAP ascertainment"],
             s7[, analysis := "Post hoc: boluses >= 1000 mL or 30 mL/kg only"],
             s8[, analysis := "Post hoc: boluses <= 500 mL only"],
             s9[, analysis := "Post hoc: clone-censor-weight for the grace window"],
             s10[, analysis := "Post hoc: phenotype-stratified propensity model"])
eff[, `:=`(outcome = fcase(grepl("^map_", term), "MAP at 60 min (mm Hg)",
                           grepl("^ne_", term),  "Norepinephrine within 6 h (risk difference)",
                           grepl("^perf_", term), "Perfusion index at 60 min (ratio)"),
           group = sub("^(map|ne|perf)_", "", term))]
eff[grepl("^perf_", term), `:=`(est = exp(est), lo = exp(lo), hi = exp(hi))]   # log ratio -> ratio
curves <- prim[grepl("^curve_", term)]
curves[, c("outcome", "modifier", "g") := tstrsplit(sub("^curve_", "", term), "_")]
curves[, value := mapply(function(v, g) grid[[v]][as.integer(g)], modifier, g)]
flow_dt <- data.table(step = vapply(flow, `[`, "", 1), n = as.integer(vapply(flow, `[`, "", 2)))

fwrite(flow_dt, file.path(cfg$out_dir, "01_flow.csv"))
fwrite(t1,      file.path(cfg$out_dir, "02_table1_by_exposure.csv"))
fwrite(t1_ph,   file.path(cfg$out_dir, "03_table_by_phenotype.csv"))
fwrite(balance, file.path(cfg$out_dir, "04_balance.csv"))
fwrite(d[bolus == 1, .N, by = .(pheno, bolus_type)][, pct := round(100 * N / sum(N), 1), by = pheno][order(pheno, -N)],
       file.path(cfg$out_dir, "15_bolus_volume_by_phenotype.csv"))
fwrite(d[, .(trial_visits = .N, map60_available_pct = round(100 * mean(!is.na(map60)), 1),
             perf60_available_pct = round(100 * mean(!is.na(perf60)), 1)), by = .(pheno, bolus)][order(pheno, bolus)],
       file.path(cfg$out_dir, "14_outcome_ascertainment_by_arm.csv"))
fwrite(eff[, .(analysis, outcome, group, est, lo, hi)], file.path(cfg$out_dir, "05_effects.csv"))

# --- revision: weight diagnostics, balance within phenotype, observed outcome counts
d[, w_tv := ipw_tv(d)]; d[, w_ph := ipw_ph(d)]
wsum <- function(w, lab) data.table(weights = lab, min = min(w), p01 = quantile(w, .01), median = median(w),
                                    p99 = quantile(w, .99), max = max(w), mean = mean(w), ess = sum(w)^2 / sum(w^2), n = length(w))
wdiag <- rbind(rbindlist(lapply(0:1, function(b) wsum(d[bolus == b, w], paste0("IPW, bolus = ", b)))),
               rbindlist(lapply(0:1, function(b) wsum(d[bolus == b, w_tv], paste0("Clone-censor-weight, bolus = ", b)))),
               rbindlist(lapply(0:1, function(b) wsum(d[bolus == b, w_ph], paste0("Phenotype-stratified IPW, bolus = ", b)))))
fwrite(wdiag, file.path(cfg$out_dir, "16_weight_diagnostics.csv"))
bal_ph <- rbindlist(lapply(levels(d$pheno), function(p) { x <- d[pheno == p]
  rbindlist(lapply(bal_vars, function(v) data.table(pheno = p, variable = v,
    smd_unweighted = smd(x[[v]], x$bolus, rep(1, nrow(x))), smd_weighted = smd(x[[v]], x$bolus, x$w)))) }))
fwrite(bal_ph, file.path(cfg$out_dir, "17_balance_by_phenotype.csv"))
bal_ph2 <- rbindlist(lapply(levels(d$pheno), function(p) { x <- d[pheno == p]
  rbindlist(lapply(bal_vars, function(v) data.table(pheno = p, variable = v,
    smd_global_ipw = smd(x[[v]], x$bolus, x$w), smd_stratified_ipw = smd(x[[v]], x$bolus, x$w_ph),
    smd_ccw = smd(x[[v]], x$bolus, x$w_tv)))) }))
fwrite(bal_ph2, file.path(cfg$out_dir, "17b_balance_by_phenotype_alt_weights.csv"))
# norepinephrine started before / during the outcome window (care after bolus initiation is part of the estimand)
d[, t_ne_rel := t_ne - t_start]
ne_t <- function(x) x[, .(trial_visits = .N,
  ne_by_15 = sum(!is.na(t_ne_rel) & t_ne_rel <= 15), ne_by_45 = sum(!is.na(t_ne_rel) & t_ne_rel <= 45),
  ne_by_75 = sum(!is.na(t_ne_rel) & t_ne_rel <= 75), ne_by_360 = sum(!is.na(t_ne_rel) & t_ne_rel <= 360),
  median_min_to_ne = as.numeric(median(t_ne_rel[!is.na(t_ne_rel) & t_ne_rel <= 360])),
  ne_by_75_with_map60 = sum(!is.na(t_ne_rel) & t_ne_rel <= 75 & !is.na(map60)), with_map60 = sum(!is.na(map60)))]
ne_tab <- rbind(d[, ne_t(.SD), by = .(pheno = as.character(pheno), bolus)], d[, ne_t(.SD), by = bolus][, pheno := "Overall"], fill = TRUE)
fwrite(ne_tab[order(pheno, bolus)], file.path(cfg$out_dir, "19_norepinephrine_timing.csv"))
nobs <- rbind(d[, .(trial_visits = .N, map60 = sum(!is.na(map60)), perf_start_and_60 = sum(perf_miss == 0 & !is.na(lperf60)),
                    ne6h = .N, visits = uniqueN(id), patients = uniqueN(pid)), by = .(pheno = as.character(pheno), bolus)],
              d[, .(pheno = "Overall", trial_visits = .N, map60 = sum(!is.na(map60)), perf_start_and_60 = sum(perf_miss == 0 & !is.na(lperf60)),
                    ne6h = .N, visits = uniqueN(id), patients = uniqueN(pid)), by = bolus])
fwrite(nobs[order(pheno, bolus)], file.path(cfg$out_dir, "18_observed_n_by_outcome.csv"))
fwrite(curves[, .(outcome, modifier, value, est, lo, hi)], file.path(cfg$out_dir, "06_effect_curves.csv"))
fwrite(d[, .(id, k, t_start, pheno, bolus, bolus_type, n_prior, ctrl_late_bolus, sbp_cur, dap, pp, map_cur, dsi,
             perf_raw, map60, perf60, ne6h, icu, w)], file.path(cfg$out_dir, "07_analytic_dataset.csv"))

f1 <- eff[analysis == "Primary (sequential trials, IPW + outcome model)" & !grepl("contrast", group)]
f1 <- f1[!is.na(est) & group %in% c("Overall", "Low pulse pressure", "Vasodilatory", "Mixed")]
f1[grepl("^Norepinephrine", outcome), `:=`(est = 100 * est, lo = 100 * lo, hi = 100 * hi)]
f1[, panel := factor(fcase(grepl("^MAP", outcome), "Difference in MAP\nat 60 min (mm Hg)",
                           grepl("^Norepinephrine", outcome), "Risk difference, norepinephrine\nwithin 6 h (percentage points)",
                           default = "Ratio of perfusion index\nat 60 min"),
                     levels = c("Difference in MAP\nat 60 min (mm Hg)", "Ratio of perfusion index\nat 60 min",
                                "Risk difference, norepinephrine\nwithin 6 h (percentage points)"))]
f1[, ref := fifelse(grepl("^Ratio", panel), 1, 0)]
f1[, group := factor(group, levels = rev(c("Overall", "Low pulse pressure", "Vasodilatory", "Mixed")))]
p1 <- ggplot(f1, aes(x = est, y = group)) +
  geom_vline(aes(xintercept = ref), linetype = "dashed", linewidth = 0.4) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4) +
  facet_wrap(~ panel, scales = "free_x", strip.position = "bottom") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.placement = "outside", panel.spacing.x = unit(1.5, "lines"))
ggsave(file.path(cfg$out_dir, "fig1_effect_by_phenotype.png"), p1, width = 10, height = 3.4, dpi = 300)
ggsave(file.path(cfg$out_dir, "fig1_effect_by_phenotype.tiff"), p1, width = 10, height = 3.4, dpi = 600, device = "tiff", compression = "lzw")
ggsave(file.path(cfg$out_dir, "fig1_effect_by_phenotype.pdf"), p1, width = 10, height = 3.4)

# Figure 6: effect of the bolus by number of boluses already received
f6 <- eff[analysis == "Primary (sequential trials, IPW + outcome model)" & grepl("^Prior boluses", group) & !is.na(est)]
f6[grepl("^Norepinephrine", outcome), `:=`(est = 100 * est, lo = 100 * lo, hi = 100 * hi)]
f6[, panel := factor(fcase(grepl("^MAP", outcome), "Difference in MAP\nat 60 min (mm Hg)",
                           grepl("^Norepinephrine", outcome), "Risk difference, norepinephrine\nwithin 6 h (percentage points)",
                           default = "Ratio of perfusion index\nat 60 min"),
                     levels = c("Difference in MAP\nat 60 min (mm Hg)", "Ratio of perfusion index\nat 60 min",
                                "Risk difference, norepinephrine\nwithin 6 h (percentage points)"))]
f6[, ref := fifelse(grepl("^Ratio", panel), 1, 0)]
f6[, group := factor(sub("Prior boluses: ", "", group), levels = rev(c("0", "1", ">=2")),
                     labels = rev(c("No prior bolus", "1 prior bolus", "2 or more prior boluses")))]
p6 <- ggplot(f6, aes(x = est, y = group)) +
  geom_vline(aes(xintercept = ref), linetype = "dashed", linewidth = 0.4) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4) +
  facet_wrap(~ panel, scales = "free_x", strip.position = "bottom") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.placement = "outside", panel.spacing.x = unit(1.5, "lines"))
ggsave(file.path(cfg$out_dir, "fig6_effect_by_prior_boluses.png"), p6, width = 10, height = 3, dpi = 300)
ggsave(file.path(cfg$out_dir, "fig6_effect_by_prior_boluses.tiff"), p6, width = 10, height = 3, dpi = 600, device = "tiff", compression = "lzw")
ggsave(file.path(cfg$out_dir, "fig6_effect_by_prior_boluses.pdf"), p6, width = 10, height = 3)

f2 <- copy(curves)
f2[outcome == "ne", `:=`(est = 100 * est, lo = 100 * lo, hi = 100 * hi)]
f2[, outcome := factor(outcome, c("map", "ne"),
                       c("Difference in MAP\nat 60 min (mm Hg)", "Risk difference, norepinephrine\n(percentage points)"))]
f2[, modifier := factor(modifier, c("dap", "pp", "dsi"),
                        c("Diastolic arterial pressure\nat trial start (mm Hg)", "Pulse pressure\nat trial start (mm Hg)",
                          "Diastolic shock index\nat trial start (bpm/mm Hg)"))]
p2 <- ggplot(f2, aes(x = value, y = est)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2) +
  geom_line(linewidth = 0.7) +
  facet_grid(outcome ~ modifier, scales = "free", switch = "both") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.placement = "outside", panel.spacing.x = unit(1.5, "lines"))
ggsave(file.path(cfg$out_dir, "fig2_effect_by_dap_pp_dsi.png"), p2, width = 11, height = 6, dpi = 300)
ggsave(file.path(cfg$out_dir, "fig2_effect_by_dap_pp_dsi.tiff"), p2, width = 11, height = 6, dpi = 600, device = "tiff", compression = "lzw")
ggsave(file.path(cfg$out_dir, "fig2_effect_by_dap_pp_dsi.pdf"), p2, width = 11, height = 6)

writeLines(c(
  "Figure 1. Effect of a crystalloid bolus started within 15 minutes of trial start on mean arterial pressure (MAP) at 60 minutes, on the perfusion index at 60 minutes (ratio, among trial-visits with perfusion index available at trial start) and on norepinephrine initiation within 6 hours, overall and by hemodynamic phenotype, in sequential emulated trials started every 15 minutes during cuff hypotension (MAP < 65 mm Hg or systolic pressure < 90 mm Hg) in emergency department patients with infection. Low pulse pressure: diastolic arterial pressure (DAP) > 50 mm Hg and pulse pressure < 40 mm Hg; vasodilatory: DAP <= 50 mm Hg and pulse pressure >= 40 mm Hg; mixed: both low. Estimates from inverse probability-weighted outcome models with g-computation, pooled across trials; bars are 95% confidence intervals from a bootstrap resampling patients.",
  "",
  "Figure 2. Effect of a crystalloid bolus across the observed range of diastolic arterial pressure, pulse pressure and diastolic shock index (heart rate divided by diastolic arterial pressure) at trial start. Lines are g-computation estimates from inverse probability-weighted models with natural-spline interaction terms, pooled across sequential emulated trials; shaded bands are 95% confidence intervals from a bootstrap resampling patients; dashed lines mark no effect."),
  file.path(cfg$out_dir, "figure_legends.txt"))

## ---- 10. Approximate Bayesian analysis -----------------------------------------
# Normal likelihood (point estimate, bootstrap SE) combined with normal priors
# (conjugate update). Priors are on the effect of the bolus (MAP, mm Hg; norepinephrine,
# risk difference) and on the between-phenotype contrasts.
step("10. Bayesian analysis")
priors <- rbind(
  data.table(outcome = "map", prior = "Skeptical",    m0 = 0, s0 = 5 / 1.96),   # 95% of mass within +/- 5 mm Hg
  data.table(outcome = "map", prior = "Weakly informative", m0 = 0, s0 = 10),
  data.table(outcome = "map", prior = "Enthusiastic", m0 = 5, s0 = 5 / 1.96),   # centered on +5 mm Hg
  data.table(outcome = "ne",  prior = "Skeptical",    m0 = 0, s0 = 0.10 / 1.96),# 95% within +/- 10 points
  data.table(outcome = "ne",  prior = "Weakly informative", m0 = 0, s0 = 0.25),
  data.table(outcome = "ne",  prior = "Enthusiastic", m0 = -0.10, s0 = 0.10 / 1.96),  # centered on 10 points fewer
  data.table(outcome = "perf", prior = "Skeptical",    m0 = 0, s0 = log(1.5) / 1.96),  # log ratio; 95% within 0.67-1.5
  data.table(outcome = "perf", prior = "Weakly informative", m0 = 0, s0 = 1),
  data.table(outcome = "perf", prior = "Enthusiastic", m0 = log(1.2), s0 = log(1.5) / 1.96))
lik <- prim[!grepl("^curve_", term)][, .(term, est, se)]
lik[, `:=`(outcome = sub("_.*", "", term), group = sub("^(map|ne|perf)_", "", term))]
post <- merge(lik, priors, by = "outcome", allow.cartesian = TRUE)
post[, `:=`(pv = 1 / (1 / s0^2 + 1 / se^2))]
post[, `:=`(pm = pv * (m0 / s0^2 + est / se^2), psd = sqrt(pv))]
post[, `:=`(cri_lo = qnorm(.025, pm, psd), cri_hi = qnorm(.975, pm, psd))]
# probabilities: MAP benefit = increase; norepinephrine benefit = decrease
post[outcome == "map", `:=`(p_any = 1 - pnorm(0, pm, psd), p_mcid = 1 - pnorm(cfg$mcid_map, pm, psd),
                            p_rope = pnorm(cfg$mcid_map, pm, psd) - pnorm(-cfg$mcid_map, pm, psd))]
post[outcome == "perf", `:=`(p_any = 1 - pnorm(0, pm, psd), p_mcid = 1 - pnorm(log(1 + cfg$mcid_perf), pm, psd),
                             p_rope = pnorm(log(1 + cfg$mcid_perf), pm, psd) - pnorm(-log(1 + cfg$mcid_perf), pm, psd))]
post[outcome == "perf", `:=`(pm = exp(pm), cri_lo = exp(cri_lo), cri_hi = exp(cri_hi), m0 = exp(m0))]   # ratio scale
post[outcome == "ne",  `:=`(p_any = pnorm(0, pm, psd), p_mcid = pnorm(-cfg$mcid_ne, pm, psd),
                            p_rope = pnorm(cfg$mcid_ne, pm, psd) - pnorm(-cfg$mcid_ne, pm, psd))]
bayes_tab <- post[, .(outcome = fcase(outcome == "map", "MAP at 60 min (mm Hg)", outcome == "ne", "Norepinephrine within 6 h (risk difference)",
                                       outcome == "perf", "Perfusion index at 60 min (ratio; prior SD on log scale)"),
                      group, prior, prior_mean = m0, prior_sd = s0, post_mean = pm, cri_lo, cri_hi,
                      p_benefit_any = p_any, p_benefit_mcid = p_mcid, p_within_rope = p_rope)]
setorder(bayes_tab, outcome, group, prior)
fwrite(bayes_tab, file.path(cfg$out_dir, "09_bayesian_posteriors.csv"))
# revision: same update for the post hoc large-bolus analysis (MAP only)
pl <- merge(s7[grepl("^map_", term), .(term, est, se, outcome = "map", group = sub("^map_", "", term))],
            priors, by = "outcome", allow.cartesian = TRUE)
pl[, pv := 1 / (1 / s0^2 + 1 / se^2)][, `:=`(pm = pv * (m0 / s0^2 + est / se^2), psd = sqrt(pv))]
pl[, `:=`(cri_lo = qnorm(.025, pm, psd), cri_hi = qnorm(.975, pm, psd), p_any = 1 - pnorm(0, pm, psd),
          p_mcid = 1 - pnorm(cfg$mcid_map, pm, psd), p_rope = pnorm(cfg$mcid_map, pm, psd) - pnorm(-cfg$mcid_map, pm, psd))]
fwrite(pl[, .(group, prior, post_mean = pm, cri_lo, cri_hi, p_benefit_any = p_any, p_benefit_mcid = p_mcid, p_within_rope = p_rope)],
       file.path(cfg$out_dir, "09b_bayesian_large_boluses_map.csv"))

# Figure 3: posterior densities of the bolus effect on MAP by phenotype (skeptical prior)
f3 <- post[outcome == "map" & group %in% c("Overall", "Low pulse pressure", "Vasodilatory", "Mixed")]
xg <- seq(-4, 5, by = 0.05)
f3d <- f3[, .(x = xg, dens = dnorm(xg, pm, psd)), by = .(group, prior)]
f3d[, group := factor(group, levels = c("Overall", "Low pulse pressure", "Vasodilatory", "Mixed"))]
f3d[, prior := factor(prior, levels = c("Skeptical", "Weakly informative", "Enthusiastic"))]
p3 <- ggplot(f3d, aes(x = x, y = dens, linetype = prior)) +
  annotate("rect", xmin = -cfg$mcid_map, xmax = cfg$mcid_map, ymin = -Inf, ymax = Inf, alpha = 0.12) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ group, ncol = 2) +
  labs(x = "Difference in MAP at 60 min (mm Hg)", y = "Posterior density", linetype = "Prior") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), legend.position = "bottom")
ggsave(file.path(cfg$out_dir, "fig3_bayesian_posteriors_map.png"), p3, width = 8, height = 6, dpi = 300)
ggsave(file.path(cfg$out_dir, "fig3_bayesian_posteriors_map.tiff"), p3, width = 8, height = 6, dpi = 600, device = "tiff", compression = "lzw")
ggsave(file.path(cfg$out_dir, "fig3_bayesian_posteriors_map.pdf"), p3, width = 8, height = 6)

cat(paste0("\nFigure 3. Posterior distributions of the effect of a crystalloid bolus on mean arterial pressure (MAP) at 60 minutes, overall and by hemodynamic phenotype, under skeptical (normal, mean 0 mm Hg, 95% of mass within +/- 5 mm Hg), weakly informative (mean 0, SD 10 mm Hg) and enthusiastic (mean +5 mm Hg, same SD as skeptical) priors, combined with a normal likelihood from the inverse probability-weighted g-computation estimate and its cluster-bootstrap standard error. The shaded band marks the region of practical equivalence (+/- ", cfg$mcid_map, " mm Hg).\n"),
    file = file.path(cfg$out_dir, "figure_legends.txt"), append = TRUE)
print(bayes_tab[prior == "Skeptical", .(outcome, group, post_mean = round(post_mean, 3),
                                        p_benefit_any = round(p_benefit_any, 3), p_benefit_mcid = round(p_benefit_mcid, 3),
                                        p_within_rope = round(p_within_rope, 3))])

## ---- 11. Phenotype transitions over time (alluvial plot) ----------------------
step("11. Phenotype transitions (alluvial)")
if (!requireNamespace("ggalluvial", quietly = TRUE)) stop("Install ggalluvial: install.packages('ggalluvial')", call. = FALSE)
tp_min <- c(0, 30, 60, 120)
tp_lab <- c("Index", "30 min", "60 min", "2 h")
st <- idx[, .(tp = tp_min), by = .(id, t_idx)][, t := t_idx + tp]
st <- merge(st, ne, by = "id", all.x = TRUE)
st <- merge(st, vis[, .(id, dep_rel)], by = "id")
n_all <- uniqueN(st$id)
st <- st[!id %in% st[tp == 0 & !is.na(t_ne) & t_ne <= t, id]]          # norepinephrine already running at index
n_ne0 <- n_all - uniqueN(st$id)
rd <- nbp[st[, .(id, tp, t)], on = "id", allow.cartesian = TRUE, nomatch = 0][tm <= t & tm > t - cfg$trial_step]
rd <- rd[order(-tm), .SD[1], by = .(id, tp)][, .(id, tp, sbp, dbp, map)]
st <- merge(st, rd, by = c("id", "tp"), all.x = TRUE)
st[, state := fcase(
  !is.na(t_ne) & t_ne <= t,                                   "Norepinephrine",
  !is.na(dep_rel) & dep_rel <= t,                             "Left ED",
  is.na(sbp),                                                 "No reading",
  map >= cfg$map_hypo & sbp >= cfg$sbp_hypo,                  "Not hypotensive",
  dbp <= cfg$dap_low & sbp - dbp >= cfg$pp_low,               "Vasodilatory",
  dbp >  cfg$dap_low & sbp - dbp <  cfg$pp_low,               "Low pulse pressure",
  dbp <= cfg$dap_low & sbp - dbp <  cfg$pp_low,               "Mixed",
  default =                                                   "Other hypotensive")]
n_pre <- uniqueN(st$id)
st <- st[!id %in% st[tp == 0 & !state %in% c("Low pulse pressure", "Vasodilatory", "Mixed"), id]]
n_oth <- n_pre - uniqueN(st$id)

# exposure for the plot: crystalloid bolus started within 30 min of index
b30 <- bol[st[tp == 0, .(id, t_idx)], on = "id", nomatch = 0][t_b >= t_idx & t_b <= t_idx + 30, unique(id)]
fin <- st[tp == 0, .(id, axis = "6 h outcome",
                     state = fifelse(!is.na(t_ne) & t_ne <= t_idx + cfg$ne_horizon, "Norepinephrine", "No norepinephrine"))]
al <- rbind(st[, .(id, axis = tp_lab[match(tp, tp_min)], state)], fin)
al[, axis := factor(axis, levels = c(tp_lab, "6 h outcome"))]
al[, exposure := factor(fifelse(id %in% b30, "Bolus within 30 min of index", "No bolus within 30 min of index"),
                        levels = c("Bolus within 30 min of index", "No bolus within 30 min of index"))]
state_lev <- c("Low pulse pressure", "Vasodilatory", "Mixed", "Other hypotensive", "Not hypotensive",
               "Norepinephrine", "No norepinephrine", "Left ED", "No reading")
al[, state := factor(state, levels = state_lev)]
al[, pct := 100 / uniqueN(id), by = exposure]           # each visit's share of its panel
pal <- c("Low pulse pressure" = "#0072B2", "Vasodilatory" = "#E69F00", "Mixed" = "#D55E00",
         "Other hypotensive" = "#CC79A7", "Not hypotensive" = "#009E73", "Norepinephrine" = "#000000",
         "No norepinephrine" = "#BBBBBB", "Left ED" = "#56B4E9", "No reading" = "#E5E5E5")

p4 <- ggplot(al, aes(x = axis, stratum = state, alluvium = id, y = pct, fill = state)) +
  ggalluvial::geom_flow(stat = "alluvium", lode.guidance = "frontback", alpha = 0.55, width = 0.3) +
  ggalluvial::geom_stratum(width = 0.3, linewidth = 0.2, color = "white") +
  scale_fill_manual(values = pal, drop = TRUE) +
  facet_wrap(~ exposure, ncol = 1) +
  labs(x = NULL, y = "Visits (%)", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), legend.position = "bottom")
ggsave(file.path(cfg$out_dir, "fig4_alluvial_phenotype_transitions.png"), p4, width = 9, height = 8, dpi = 300)
ggsave(file.path(cfg$out_dir, "fig4_alluvial_phenotype_transitions.tiff"), p4, width = 9, height = 8, dpi = 600, device = "tiff", compression = "lzw")
ggsave(file.path(cfg$out_dir, "fig4_alluvial_phenotype_transitions.pdf"), p4, width = 9, height = 8)

# phenotype stability: transitions from index to 30 min, and agreement among visits still hypotensive
w0 <- dcast(st[tp %in% c(0, 30)], id ~ tp, value.var = "state"); setnames(w0, c("id", "s0", "s30"))
trans <- w0[, .N, by = .(index = s0, at_30_min = s30)][order(index, -N)]
trans[, pct_of_index := round(100 * N / sum(N), 1), by = index]
still <- w0[s30 %in% c("Low pulse pressure", "Vasodilatory", "Mixed")]
stab  <- data.table(measure = c("Visits with index hypotension", "  excluded: norepinephrine running at index",
                                "  excluded: index reading not classifiable", "Visits in alluvial cohort",
                                "Still hypotensive with a phenotype at 30 min", "  same phenotype at 30 min (%)",
                                "No cuff reading at 30 min", "Not hypotensive at 30 min"),
                    value = c(n_all, n_ne0, n_oth, uniqueN(st$id), nrow(still), round(100 * mean(still$s0 == still$s30), 1),
                              w0[s30 == "No reading", .N], w0[s30 == "Not hypotensive", .N]))
fwrite(trans, file.path(cfg$out_dir, "10_phenotype_transitions_index_to_30min.csv"))
fwrite(stab,  file.path(cfg$out_dir, "11_phenotype_stability.csv"))
cat("\nFigure 4. Transitions of hemodynamic phenotype from the first hypotensive cuff reading (index) to 30 minutes, 60 minutes and 2 hours, and norepinephrine initiation within 6 hours, in visits that did and did not receive a crystalloid bolus within 30 minutes of index. Each band is a group of visits; band width is proportional to the number of visits. At each time point the state is taken from the latest cuff reading in the previous 15 minutes; visits on norepinephrine or discharged from the emergency department are shown as such from that time onward.\n",
    file = file.path(cfg$out_dir, "figure_legends.txt"), append = TRUE)
print(stab)

## ---- 12. MAP and perfusion-index trajectories around trial start ------------
step("12. Trajectories (weighted, cluster bootstrap)")
bins <- seq(-30, 120, by = 15)
get_traj <- function(src) {
  j <- src[d[, .(id, k, t_start)], on = "id", allow.cartesian = TRUE, nomatch = 0][, rel := t - t_start]
  j <- j[rel >= min(bins) & rel < max(bins) + 15][, bin := floor(rel / 15) * 15 + 7.5]   # bin midpoints
  j[, .(value = median(value)), by = .(id, k, bin)]
}
traj <- rbind(get_traj(nbp[, .(id, t = tm, value = map)])[, outcome := "map"],
              get_traj(nm[var == "perf", .(id, t = t_rel, value = log(value))])[, outcome := "lperf"])
traj_means <- function(x) {             # x: trial-visits with weights; may contain duplicated (id, k) rows
  x <- x[, .(id, k, bolus, pheno = as.character(pheno), w)]
  j <- merge(traj, x, by = c("id", "k"), allow.cartesian = TRUE)
  rbind(j[, .(m = weighted.mean(value, w)), by = .(outcome, bin, bolus, pheno)],
        j[, .(m = weighted.mean(value, w), pheno = "Overall"), by = .(outcome, bin, bolus)])
}
d[, w := ipw(d)]
tm_pt <- traj_means(d)
tm_bs <- rbindlist(lapply(seq_len(cfg$n_boot_traj), function(b) {
  xb <- cluster_sample(d)
  xb[, w := tryCatch(ipw(xb), error = function(e) NA_real_)]
  if (anyNA(xb$w)) return(NULL)
  traj_means(xb)[, rep := b]
}))
tm_ci <- tm_bs[, .(lo = quantile(m, .025), hi = quantile(m, .975)), by = .(outcome, bin, bolus, pheno)]
tm_all <- merge(tm_pt, tm_ci, by = c("outcome", "bin", "bolus", "pheno"))
# difference bolus - no bolus per bin
dif_pt <- dcast(tm_pt, outcome + bin + pheno ~ bolus, value.var = "m")[, diff := `1` - `0`]
dif_bs <- dcast(tm_bs, outcome + bin + pheno + rep ~ bolus, value.var = "m")[, diff := `1` - `0`]
dif <- merge(dif_pt[, .(outcome, bin, pheno, diff)],
             dif_bs[, .(lo = quantile(diff, .025, na.rm = TRUE), hi = quantile(diff, .975, na.rm = TRUE)),
                    by = .(outcome, bin, pheno)], by = c("outcome", "bin", "pheno"))
tm_all[outcome == "lperf", `:=`(m = exp(m), lo = exp(lo), hi = exp(hi))]
dif[outcome == "lperf", `:=`(diff = exp(diff), lo = exp(lo), hi = exp(hi))]
fwrite(tm_all, file.path(cfg$out_dir, "12_trajectories_by_arm.csv"))
fwrite(dif,    file.path(cfg$out_dir, "13_trajectory_differences.csv"))   # perf: ratio of geometric means

f5 <- copy(tm_all)
f5[, arm := factor(fifelse(bolus == 1, "Bolus", "No bolus"), levels = c("Bolus", "No bolus"))]
f5[, pheno := factor(pheno, levels = c("Overall", "Low pulse pressure", "Vasodilatory", "Mixed"))]
f5[, outcome := factor(outcome, c("map", "lperf"), c("MAP (mm Hg)", "Perfusion index\n(geometric mean)"))]
p5 <- ggplot(f5, aes(x = bin, y = m, color = arm, fill = arm)) +
  annotate("rect", xmin = 0, xmax = cfg$grace, ymin = -Inf, ymax = Inf, alpha = 0.08) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1) +
  scale_color_manual(values = c("Bolus" = "#0072B2", "No bolus" = "#7F7F7F")) +
  scale_fill_manual(values = c("Bolus" = "#0072B2", "No bolus" = "#7F7F7F")) +
  scale_x_continuous(breaks = seq(0, 120, 60)) +
  facet_grid(outcome ~ pheno, scales = "free_y", switch = "y") +
  labs(x = "Minutes from trial start", y = NULL, color = NULL, fill = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.placement = "outside", legend.position = "bottom",
        panel.spacing.x = unit(1.5, "lines"))
ggsave(file.path(cfg$out_dir, "fig5_trajectories_map_perf.png"), p5, width = 11, height = 6, dpi = 300)
ggsave(file.path(cfg$out_dir, "fig5_trajectories_map_perf.tiff"), p5, width = 11, height = 6, dpi = 600, device = "tiff", compression = "lzw")
ggsave(file.path(cfg$out_dir, "fig5_trajectories_map_perf.pdf"), p5, width = 11, height = 6)
cat("\nFigure 5. Inverse probability-weighted trajectories of mean arterial pressure (MAP) and perfusion index (geometric mean) from 30 minutes before to 2 hours after trial start in trial-visits with and without a crystalloid bolus started within 15 minutes of trial start (shaded column), overall and by hemodynamic phenotype at trial start. Values are medians within 15-minute bins plotted at bin midpoints; bands are 95% confidence intervals from a bootstrap resampling patients.\n",
    file = file.path(cfg$out_dir, "figure_legends.txt"), append = TRUE)

cat("\nFigure 6. Effect of a crystalloid bolus started within 15 minutes of trial start on mean arterial pressure (MAP) at 60 minutes, perfusion index at 60 minutes (ratio) and norepinephrine initiation within 6 hours, by number of crystalloid boluses already received before trial start (pre-specified effect modifier). Estimates from inverse probability-weighted outcome models with a bolus-by-prior-boluses interaction and g-computation, pooled across sequential emulated trials; bars are 95% confidence intervals from a bootstrap resampling patients.\n",
    file = file.path(cfg$out_dir, "figure_legends.txt"), append = TRUE)
## ---- 13. Extra figures A and B (from objects in memory) ---------------------
step("13. Extra figures A and B")
ph_lev  <- c("Overall", "Low pulse pressure", "Vasodilatory", "Mixed")
pal_ph  <- c("Overall" = "#000000", "Low pulse pressure" = "#0072B2", "Vasodilatory" = "#E69F00", "Mixed" = "#D55E00")
arm_col <- c("Bolus" = "#0072B2", "No bolus" = "#7F7F7F")

# Figure A: change in MAP at 60 min vs MAP at trial start, by arm and phenotype (weighted smooths)
da <- d[!is.na(map60)]
da[, `:=`(dmap = map60 - map_cur,
          arm = factor(fifelse(bolus == 1, "Bolus", "No bolus"), levels = c("Bolus", "No bolus")))]
fa <- rbind(da[, .(pheno = as.character(pheno), map_cur, dmap, arm, w)],
            da[, .(pheno = "Overall", map_cur, dmap, arm, w)])
fa[, pheno := factor(pheno, levels = ph_lev)]
fa <- fa[map_cur >= quantile(map_cur, .02) & map_cur <= quantile(map_cur, .98)]   # trim sparse tails
pA <- ggplot(fa, aes(x = map_cur, y = dmap, color = arm, fill = arm)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(alpha = 0.08, size = 0.6, position = position_jitter(width = 0.3, height = 0)) +
  geom_smooth(aes(weight = w), method = "loess", span = 0.9, linewidth = 0.9, alpha = 0.2) +
  scale_color_manual(values = arm_col) + scale_fill_manual(values = arm_col) +
  coord_cartesian(ylim = c(-15, 35)) +
  facet_wrap(~ pheno, nrow = 1) +
  labs(x = "MAP at trial start (mm Hg)", y = "Change in MAP at 60 min (mm Hg)", color = NULL, fill = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), legend.position = "bottom", panel.spacing.x = unit(1.2, "lines"))
ggsave(file.path(cfg$out_dir, "figA_map_response_by_baseline.png"), pA, width = 11, height = 4, dpi = 300)
ggsave(file.path(cfg$out_dir, "figA_map_response_by_baseline.tiff"), pA, width = 11, height = 4, dpi = 600, device = "tiff", compression = "lzw")
ggsave(file.path(cfg$out_dir, "figA_map_response_by_baseline.pdf"), pA, width = 11, height = 4)

# Figure B: posterior probability of benefit above each threshold (skeptical prior), MAP and perfusion index
pb <- post[prior == "Skeptical" & group %in% ph_lev & outcome %in% c("map", "perf")]
pb[outcome == "perf", pm := log(pm)]                       # back to the log scale (post stores ratios)
grid_b <- rbind(data.table(outcome = "map",  t = seq(0, 5, by = 0.05)),
                data.table(outcome = "perf", t = seq(0, 30, by = 0.25)))   # perf threshold: % increase
fb <- merge(pb[, .(outcome, group, pm, psd)], grid_b, by = "outcome", allow.cartesian = TRUE)
fb[outcome == "map",  p := 1 - pnorm(t, pm, psd)]
fb[outcome == "perf", p := 1 - pnorm(log(1 + t / 100), pm, psd)]
fb[, group := factor(group, levels = ph_lev)]
fb[, panel := factor(outcome, c("map", "perf"),
                     c("Increase in MAP at 60 min greater than (mm Hg)",
                       "Increase in perfusion index at 60 min greater than (%)"))]
mcid_lines <- data.table(panel = factor(levels(fb$panel), levels = levels(fb$panel)), x = c(cfg$mcid_map, 100 * cfg$mcid_perf))
pB <- ggplot(fb, aes(x = t, y = p, color = group)) +
  geom_vline(data = mcid_lines, aes(xintercept = x), linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = pal_ph) +
  scale_y_continuous(limits = c(0, 1), labels = function(z) paste0(100 * z, "%")) +
  facet_wrap(~ panel, scales = "free_x", strip.position = "bottom") +
  labs(x = NULL, y = "Posterior probability", color = NULL) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.placement = "outside",
        legend.position = "bottom", panel.spacing.x = unit(1.5, "lines"))
ggsave(file.path(cfg$out_dir, "figB_probability_of_benefit.png"), pB, width = 10, height = 4.2, dpi = 300)
ggsave(file.path(cfg$out_dir, "figB_probability_of_benefit.tiff"), pB, width = 10, height = 4.2, dpi = 600, device = "tiff", compression = "lzw")
ggsave(file.path(cfg$out_dir, "figB_probability_of_benefit.pdf"), pB, width = 10, height = 4.2)

cat("\nFigure A. Change in mean arterial pressure (MAP) from trial start to 60 minutes as a function of MAP at trial start, in trial-visits with and without a crystalloid bolus started within 15 minutes, overall and by hemodynamic phenotype. Points are individual trial-visits; lines are inverse probability-weighted locally estimated smooths with 95% bands (descriptive; bands do not account for repeated trials per patient). The 2nd and 98th percentiles of MAP at trial start bound the x axis.\n",
    "\nFigure B. Posterior probability that a crystalloid bolus increases mean arterial pressure (MAP) at 60 minutes by more than each threshold (left) and the perfusion index at 60 minutes by more than each relative threshold (right), overall and by hemodynamic phenotype, under a skeptical prior. Dashed lines mark the pre-specified minimal clinically important differences (2 mm Hg and 10%).\n",
    file = file.path(cfg$out_dir, "figure_legends.txt"), append = TRUE)

print(flow_dt)
print(eff[grepl("^Primary", analysis), .(outcome, group, est = round(est, 3), lo = round(lo, 3), hi = round(hi, 3))])
message("Done. Outputs in: ", normalizePath(cfg$out_dir))

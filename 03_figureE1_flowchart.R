# =====================================================================
# MC-MED — Figure E1: flow of visits and emulated trials
# Reads 01_flow.csv and 18_observed_n_by_outcome.csv from the output folder
# (base R + grid only). No title inside the figure; legend is separate.
# =====================================================================
library(grid)
home    <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else path.expand("~")
out_dir <- Sys.getenv("MCMED_OUT_DIR", unset = file.path(home, "Desktop", "mimiciv", "analysis_dap_pp_bolus_seqtrials"))
if (!dir.exists(out_dir)) out_dir <- getwd()               # fallback: files in the working directory
n_mcmed <- 118385                                           # adult visits in MC-MED v1.0.1

fl <- read.csv(file.path(out_dir, "01_flow.csv"), stringsAsFactors = FALSE)
ob <- read.csv(file.path(out_dir, "18_observed_n_by_outcome.csv"), stringsAsFactors = FALSE)
g  <- function(pat) fl$n[grep(pat, fl$step, fixed = TRUE)][1]
f  <- function(x) formatC(x, big.mark = ",", format = "d")

n_inf <- g("Infection visits");               n_hyp <- g("Visits with cuff hypotension")
n_cand <- g("Candidate trial-visits");         n_ne  <- g("no norepinephrine started")
n_ed  <- g("still in ED");                     n_cuff <- g("cuff reading in previous")
n_hypo <- g("currently hypotensive");          n_wash <- g("no bolus in previous 30 min")
n_oth <- g("phenotype other");                 n_an  <- g("Analytic trial-visits")
n_vis <- g("distinct visits");                 n_pat <- g("distinct patients")
n_b   <- g("trial-visits with bolus");         n_bv  <- g("distinct visits ever exposed")
o1 <- ob[ob$pheno == "Overall" & ob$bolus == 1, ]; o0 <- ob[ob$pheno == "Overall" & ob$bolus == 0, ]
ph <- sapply(c("Low pulse pressure", "Vasodilatory", "Mixed"), function(p) sum(ob$trial_visits[ob$pheno == p]))

box <- function(x, y, w, h, lab, bold_first = TRUE, just = "centre") {
  grid.rect(unit(x, "npc"), unit(y, "npc"), unit(w, "npc"), unit(h, "npc"), gp = gpar(fill = "white", col = "black", lwd = 0.8))
  grid.text(lab, x = if (just == "left") unit(x - w / 2, "npc") + unit(2.5, "mm") else unit(x, "npc"), y = unit(y, "npc"),
            just = just, gp = gpar(fontsize = 8.5, lineheight = 1.15))
}
arr <- function(x0, y0, x1, y1) grid.lines(unit(c(x0, x1), "npc"), unit(c(y0, y1), "npc"),
                                           arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
                                           gp = gpar(fill = "black", lwd = 0.8))
draw <- function() {
  grid.newpage()
  xm <- 0.33; xs <- 0.79; wm <- 0.46; ws <- 0.40
  ys <- c(0.95, 0.83, 0.71, 0.585, 0.355)
  box(xm, ys[1], wm, 0.055, paste0("Adult ED visits in MC-MED\nn = ", f(n_mcmed)))
  box(xm, ys[2], wm, 0.065, paste0("Visits with infection\n(ICD-10 sepsis code, or blood culture and antibiotic\nwithin 24 h)  n = ", f(n_inf)))
  box(xm, ys[3], wm, 0.065, paste0("Visits with cuff hypotension (MAP < 65 or SBP < 90 mm Hg)\nwithin 6 h of rooming\nn = ", f(n_hyp)))
  box(xm, ys[4], wm, 0.075, paste0("Candidate trial-visits\n(a trial every 15 min from the first hypotensive\nreading up to 3 h later)  n = ", f(n_cand)))
  box(xm, ys[5], wm, 0.11, paste0("Eligible trial-visits  n = ", f(n_an), "\n", f(n_vis), " visits, ", f(n_pat), " patients\n",
                                  "Low pulse pressure ", f(ph[1]), "; vasodilatory ", f(ph[2]), "; mixed ", f(ph[3])))
  for (i in 1:4) arr(xm, ys[i] - c(0.0275, 0.0325, 0.0325, 0.0375)[i], xm, ys[i + 1] + c(0.0325, 0.0325, 0.0375, 0.055)[i])
  # exclusion boxes
  ex <- list(c(0.89, paste0("No infection  n = ", f(n_mcmed - n_inf))),
             c(0.77, paste0("No cuff hypotension within 6 h  n = ", f(n_inf - n_hyp))))
  for (e in ex) { y <- as.numeric(e[1]); box(xs, y, ws, 0.035, e[2]); arr(xm, y, xs - ws / 2, y) }
  ytr <- 0.49
  box(xs, ytr, ws, 0.15, just = "left", paste0(
    "Trial-visits excluded  n = ", f(n_cand - n_an), "\n",
    "  Norepinephrine started before trial start: ", f(n_cand - n_ne), "\n",
    "  Left the ED before trial start: ", f(n_ne - n_ed), "\n",
    "  No cuff reading in the previous 15 min: ", f(n_ed - n_cuff), "\n",
    "  Latest reading not hypotensive: ", f(n_cuff - n_hypo), "\n",
    "  Crystalloid bolus in the previous 30 min: ", f(n_hypo - n_wash), "\n",
    "  DAP > 50 and PP \u2265 40 mm Hg: ", f(n_oth)))
  arr(xm, ytr, xs - ws / 2, ytr)
  # arms
  ya <- 0.15; xa <- c(0.27, 0.73); wa <- 0.4
  grid.lines(unit(c(xm, xm), "npc"), unit(c(ys[5] - 0.055, 0.255), "npc"), gp = gpar(lwd = 0.8))
  grid.lines(unit(xa, "npc"), unit(c(0.255, 0.255), "npc"), gp = gpar(lwd = 0.8))
  for (i in 1:2) arr(xa[i], 0.255, xa[i], ya + 0.075)
  box(xa[1], ya, wa, 0.15, paste0("Bolus started within 15 min\n", f(n_b), " trial-visits (", f(n_bv), " visits)\n\n",
                                  "MAP at 60 min available: ", f(o1$map60), "\nPerfusion index at start and\n60 min available: ", f(o1$perf_start_and_60)))
  box(xa[2], ya, wa, 0.15, paste0("No bolus within 15 min\n", f(n_an - n_b), " trial-visits (", f(o0$visits), " visits)\n\n",
                                  "MAP at 60 min available: ", f(o0$map60), "\nPerfusion index at start and\n60 min available: ", f(o0$perf_start_and_60)))
}
png(file.path(out_dir, "figS1_flowchart.png"), width = 7.5, height = 9, units = "in", res = 300); draw(); invisible(dev.off())
tiff(file.path(out_dir, "figS1_flowchart.tiff"), width = 7.5, height = 9, units = "in", res = 600, compression = "lzw"); draw(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "figS1_flowchart.pdf"), width = 7.5, height = 9); draw(); invisible(dev.off())
cat("\nFigure E1. Flow of visits and emulated trials. A new trial started at the first hypotensive cuff reading and every 15 minutes thereafter, up to 3 hours later; eligibility was assessed with information available at trial start, and exclusions are applied in the order listed. A visit can contribute to several trials and a patient to several visits, so visits and patients can appear in both arms. DAP, diastolic arterial pressure; ED, emergency department; MAP, mean arterial pressure; PP, pulse pressure; SBP, systolic blood pressure.\n",
    file = file.path(out_dir, "figure_legends.txt"), append = TRUE)
message("Figure E1 saved in ", normalizePath(out_dir))

# scripts/sims/sim_canonicals_run_all.R
# Run comprehensive simulations for all three DGPs

source(here::here("scripts/sims/canonical_run_comprehensive_simulation.R"))

cat("\n", rep("=", 80), "\n", sep = "")
cat("RUNNING ALL COMPREHENSIVE SIMULATIONS\n")
cat(rep("=", 80), "\n\n", sep = "")

# ACIC Simulation
cat(">>> Starting ACIC comprehensive simulation...\n")
res_acic <- run_comprehensive_sim(
  sim_type = "acic",
  n_iter = 100,
  seed_start = 1000
)

# Hainmueller Simulation
cat("\n>>> Starting Hainmueller comprehensive simulation...\n")
res_hain <- run_comprehensive_sim(
  sim_type = "hainmueller",
  n_iter = 100,
  seed_start = 1000
)

# Kang & Schafer Simulation
cat("\n>>> Starting Kang & Schafer comprehensive simulation...\n")
res_kang <- run_comprehensive_sim(
  sim_type = "kang",
  n_iter = 100,
  seed_start = 1000
)

cat("\n", rep("=", 80), "\n", sep = "")
cat("ALL SIMULATIONS COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")

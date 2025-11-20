# scripts/sims/run_local_all.R
# Local runner: iterate run_single_iteration.R for chosen sims.
# Usage:
#   Rscript scripts/sims/run_local_all.R [n_iter] [methods...]
# Examples:
#   Rscript scripts/sims/run_local_all.R                 # 5 iters, all methods
#   Rscript scripts/sims/run_local_all.R 3 diff twang    # 3 iters, only diff & twang
#   Rscript scripts/sims/run_local_all.R 10 "ps,dr"      # 10 iters, PS + doubly-robust groups

args <- commandArgs(trailingOnly = TRUE)

# default iterations
n_iter <- 5L
if (length(args) >= 1L && nzchar(args[1])) {
  maybe_int <- suppressWarnings(as.integer(args[1]))
  if (!is.na(maybe_int)) {
    n_iter <- maybe_int
    args <- args[-1]
  }
}

# parse optional methods (space- and/or comma-separated)
parse_methods <- function(x) {
  if (length(x) == 0) return(character(0))
  m <- unlist(strsplit(x, split = ","))
  m <- trimws(tolower(m))
  m[nzchar(m)]
}
methods <- parse_methods(args)

method_suffix <- if (length(methods) > 0) {
  paste("", paste(shQuote(methods), collapse = " "), sep = " ")
} else {
  ""
}

run_many <- function(sim, n_iter = 5L, method_suffix = "") {
  message(sprintf(">>> %s: running %d iterations%s",
                  sim, n_iter,
                  if (nzchar(method_suffix)) paste0(" with methods ", paste(methods, collapse = ", ")) else " (all methods)"))
  for (i in seq_len(n_iter)) {
    cmd <- sprintf("Rscript scripts/sims/run_single_iteration.R %s %d%s", sim, i, method_suffix)
    status <- system(cmd)
    if (status != 0) stop(sprintf("Failed on %s iteration %d", sim, i))
  }
}

sims <- c("acic","hainmueller","kang","toy")
for (s in sims) run_many(s, n_iter = n_iter, method_suffix = method_suffix)

message("All local runs complete.")

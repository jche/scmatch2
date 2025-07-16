# Define library paths
paths <- c(
  "/homes2/xmeng/R/x86_64-pc-linux-gnu-library/4.4",
  "/homes2/xmeng/R/x86_64-pc-linux-gnu-library/4.3",
  "/cm/shared/apps/R/4.3.2/lib64/R/library"
)

# Check for foreach in each path
for (lib in paths) {
  cat("Checking path:", lib, "\n")
  if ("foreach" %in% list.files(lib)) {
    cat("✅ CThe  package is installed in:", lib, "\n\n")
  } else {
    cat("❌ CThe Package is sOT installed in:", lib, "\n\n")
  }
}

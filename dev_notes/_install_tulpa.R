# Rebuild + install upstream tulpa after the ABI bump (21 -> 22) and the
# slope-only RE engine change. tulpaObs LinkingTo: tulpa picks up the new
# headers only after tulpa is installed (not merely load_all'd).
#
# CLEAN rebuild: header-only changes (model_data.h, log_post_generic_impl.h)
# do not reliably trigger recompilation of the .cpp/.o that include them, so
# wipe the compiled objects first or the ABI getter + AD gradient stay stale.
pkgbuild::clean_dll("C:/Users/Gilles Colling/Documents/dev/tulpa")
devtools::install(
  "C:/Users/Gilles Colling/Documents/dev/tulpa",
  quick = TRUE, upgrade = FALSE, reload = FALSE
)
cat("=== tulpa installed: ABI",
    tryCatch(tulpa:::.tulpa_abi_version(), error = function(e) "n/a"), "===\n")

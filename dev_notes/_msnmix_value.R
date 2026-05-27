db <- tools::Rd_db("spAbundance")
rd <- db[["msNMix.Rd"]]
tmp <- tempfile(fileext = ".txt")
tools::Rd2txt(rd, out = tmp)
txt <- readLines(tmp)
# Print from "Value" onward and the Examples
val <- grep("^_V_a_l_u_e", txt)
ex  <- grep("^_E_x_a_m_p_l_e", txt)
cat("=== VALUE ===\n")
if (length(val)) cat(paste(txt[val:(min(ex, length(txt)) - 1)], collapse = "\n"))
cat("\n\n=== EXAMPLES ===\n")
if (length(ex)) cat(paste(txt[ex:length(txt)], collapse = "\n"))

.onAttach <- function(libname, pkgname) {
  version <- tryCatch(
    as.character(utils::packageVersion(pkgname)),
    error = function(e) "unknown"
  )

  packageStartupMessage(
    "CLIPER : Version ", version, "\n",
    "CLusterIng-enabled Peak-to-gEne Regression for peak-to-gene fine-mapping\n"
  )
}

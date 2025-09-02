# Declare global variables here to avoid warnings during package testing
utils::globalVariables(c(
  ".SD", # Special object of data.table
  ".N", # Special object of data.table
  ":=", # Special object of data.table
  ".", # Special object of data.table
  "Date", # Variable created dynamically in create_data
  "Cod", # Variable created dynamically in IDW
  "X", # Variable created dynamically in IDW
  "Y", # Variable created dynamically in IDW
  "Z", # Variable created dynamically in IDW
  "ID", # Variable created dynamically in IDW
  "sum_d", # Variable created dynamically in IDW
  "sum_n", # Variable created dynamically in IDW
  "value", # Variable created dynamically in IDW
  "x", # Variable created dynamically in IDW
  "y", # Variable created dynamically in IDW
  "observed", # Variable created dynamically in Intr_Connections_module
  "estimated", # Variable created dynamically in Intr_Connections_module
  "Obs", # Variable created dynamically in Intr_Connections_module
  "Sim", # Variable created dynamically in Intr_Connections_module,
  "var", # Variable created dynamically in IDW and Cressman
  "..features_ff", # Variable created dynamically in RFmerge and RFplus
  "residuals", # Variable created dynamically in RFplus
  "sim" # Variable created dynamically in RFplus
))

.onAttach <- function(libname, pkgname) {
  installed_version <- as.character(utils::packageVersion(pkgname))

  cran_version <- tryCatch(
    {
      available <- utils::available.packages()
      if (pkgname %in% rownames(available)) {
        available[pkgname, "Version"]
      } else {
        NA
      }
    },
    error = function(e) NA
  )

  # Mostrar mensajes
  packageStartupMessage("+--------------------------------------------+")
  packageStartupMessage(paste0(
    "| Installed version of ",
    pkgname,
    ": ",
    installed_version
  ))

  if (!is.na(cran_version)) {
    packageStartupMessage(paste0("| Version available on CRAN: ", cran_version))

    if (utils::compareVersion(installed_version, cran_version) < 0) {
      packageStartupMessage("|")
      packageStartupMessage("| *** A new version is available ***")
      packageStartupMessage("| To update, run:")
      packageStartupMessage(paste0("| install.packages(\"", pkgname, "\")"))
    }
  } else {
    packageStartupMessage(
      "| The version available on CRAN could not be consulted."
    )
  }

  packageStartupMessage("+--------------------------------------------+")
}

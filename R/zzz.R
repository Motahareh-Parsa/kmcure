.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("The 'kmcure' R package (version = ", packageVersion('kmcure'), ") is loaded :) \nYou can run its example by using 'example(kmcure)' command!"))
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage("********************************************************************************")
  packageStartupMessage("                       Welcome to the package 'fEGarch'!")
  packageStartupMessage("********************************************************************************")
  packageStartupMessage("")
  packageStartupMessage("Please report any possible errors and bugs to dominik.schulz@uni-paderborn.de.")
  packageStartupMessage("")
  packageStartupMessage("If you are using this package for your publication, please consider citing")
  packageStartupMessage('this software as indicated by the command citation("fEGarch").')
  packageStartupMessage("")
  packageStartupMessage("Thank you.")
  packageStartupMessage("")
  packageStartupMessage("********************************************************************************")
}

.onUnload <- function(libpath) {
  library.dynam.unload("fEGarch", libpath)
}

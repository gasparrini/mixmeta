###
### R routines for the R package mixmeta (c)
#
.onAttach <-
function(lib, pkg) {
#
################################################################################
#
  meta <- packageDescription("mixmeta")
  attachmsg <- paste("This is mixmeta ",meta$Version,
    ". For an overview type: help('mixmeta-package').",sep="")
  packageStartupMessage(attachmsg, domain = NULL, appendLF = TRUE)
}

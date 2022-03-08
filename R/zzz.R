#' @imclude utils.R

.onLoad <- function(libname, pkgname) {
  if (! requireNamespace("rjd3sts", quietly=T)) stop("Loading rjd3 libraries failed")

  result <- .jpackage(pkgname, lib.loc=libname)
  if (!result) stop("Loading java packages failed")

  proto.dir <- system.file("proto", package = pkgname)
  RProtoBuf::readProtoFiles2(protoPath = proto.dir)
  
  # reload extractors
  .jcall("demetra/information/InformationExtractors", "V", "reloadExtractors")
}



.onLoad <- function(libname, pkgname) {
  suppressMessages(require(rjd3sts, quietly = T))
  result <- .jpackage(pkgname, lib.loc=libname)
  if (!result) stop("Loading java packages failed")

  proto.dir <- system.file("proto", package = pkgname)
  readProtoFiles2(protoPath = proto.dir)
  
  # reload extractors
  .jcall("demetra/information/InformationExtractors", "V", "reloadExtractors")
}


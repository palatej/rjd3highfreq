#' @imclude utils.R
NULL

#' Perform an STL like (based on Loess) decomposition on any periodicity
#'
#' @param y input time series.
#' @param period period, any positive real number.
#' @param multiplicative Boolean indicating if the decomposition mode is multiplicative (TRUE).
#' @param swindow length of seasonal filter.
#' @param twindow length of trend filter.
#' @param robust Boolean: in outer loop robust weights for irregular.
#'
#' @return
#' @export
#'
#' @examples
stl<-function(y, period, multiplicative=TRUE, swindow=7, twindow=0, robust=TRUE){
  jrslt<-.jcall("demetra/stl/r/StlDecomposition", "Ldemetra/math/matrices/Matrix;", "process", as.numeric(y), as.integer(period), multiplicative, as.integer(swindow), as.integer(twindow), robust)
  m<-rjd3toolkit::matrix_jd2r(jrslt)
  decomposition<-list(
    y=m[,1],
    sa=m[,2],
    t=m[,3],
    s=m[,4],
    i=m[,5]
  )
  parameters<-list(
    multiplicative=multiplicative, 
    swindow=swindow, 
    twindow=twindow, 
    robust=robust
  )

  return(structure(list(
    decomposition=decomposition,
    parameters=parameters),
    class="JDSTL"))
}

#' Fit a Loess regression.
#'
#' @param y input time series.
#' @param window 
#' @param degree 
#' @param jump 
#'
#' @return
#' @export
#'
#' @examples
loess<-function(y, window, degree=1, jump=1){
  if (degree != 0 && degree != 1)
    stop("Unsupported degree")
  if (jump <1)
    stop("jump should be greater then 0")
  return (.jcall("demetra/r/StlDecomposition", "[D", "loess", as.numeric(y), as.integer(window), as.integer(degree), as.integer(jump)))
}


#' @imclude utils.R
NULL

ucm_extract<-function(jrslt, cmp){
  path<-paste0("ucarima.component(", cmp,")")
  return (arima_extract(jrslt, path))
}

arima_extract<-function(jrslt, path){
  str<-rjd3toolkit:::proc_str(jrslt, paste0(path, ".name"))
  ar<-rjd3toolkit:::proc_vector(jrslt, paste0(path, ".ar"))
  delta<-rjd3toolkit:::proc_vector(jrslt, paste0(path, ".delta"))
  ma<-rjd3toolkit:::proc_vector(jrslt, paste0(path, ".ma"))
  var<-rjd3toolkit:::proc_numeric(jrslt, paste0(path, ".var"))
  return (rjd3modelling::arima.model(str, ar,delta,ma,var))
}



#' Title
#'
#' @param y 
#' @param period Period of the seasonal component
#' @param adjust True if an actual fractional airline model is used. False if the period is rounded to the nearest integer
#' @param sn Signal/noise decomposition. The signal is the seasonally adjusted series and the noise the seasonal component
#' @param stde True if standard deviations of the components must be computed. In some cases (memory limits), it is currently not possible to compute them
#' @param nbcasts Number of backcasts
#' @param nfcasts Number of forecasts
#'
#' @return
#' @export
#'
#' @examples
fractionalAirlineDecomposition<-function(y, period, sn=F, stde=F, nbcasts=0, nfcasts=0){
  checkmate::assertNumeric(y, null.ok = F)
  checkmate::assertNumeric(period, len = 1, null.ok = F)
  checkmate::assertLogical(sn, len = 1, null.ok = F)
  jrslt<-.jcall("demetra/highfreq/r/FractionalAirlineProcessor", "Ljdplus/highfreq/FractionalAirlineDecomposition;", "decompose", as.numeric(y), 
                period, sn, stde, as.integer(nbcasts), as.integer(nfcasts))
  return (jd2r_fractionalAirlineDecomposition(jrslt, sn, stde))  
}

#' Title
#'
#' @param y 
#' @param periods 
#' @param ndiff 
#' @param stde 
#' @param nbcasts 
#' @param nfcasts 
#'
#' @return
#' @export
#'
#' @examples
multiAirlineDecomposition<-function(y, periods, ndiff=2, ar=F, stde=F, nbcasts=0, nfcasts=0){
  if (length(periods) == 1){
    return (fractionalAirlineDecomposition(y, periods, stde=stde, nbcasts = nbcasts, nfcasts = nfcasts))
  }
  checkmate::assertNumeric(y, null.ok = F)
  
  jrslt<-.jcall("demetra/highfreq/r/FractionalAirlineProcessor", "Ljdplus/highfreq/FractionalAirlineDecomposition;", "decompose", as.numeric(y), 
                .jarray(periods), as.integer(ndiff), ar, stde, as.integer(nbcasts), as.integer(nfcasts))
  
  if (length(periods) == 1){
    return (jd2r_fractionalAirlineDecomposition(jrslt, F, stde))
  }else{
    return (jd2r_multiAirlineDecomposition(jrslt, stde))
  }
}

#' Title
#'
#' @param y 
#' @param periods 
#' @param x 
#' @param mean 
#' @param outliers 
#' @param criticalValue 
#'
#' @return
#' @export
#'
#' @examples
fractionalAirlineEstimation<-function(y, periods, x = NULL, ndiff=2, ar=F, mean = FALSE, outliers=NULL, criticalValue=6, precision=1e-12, approximateHessian=F){
  checkmate::assertNumeric(y, null.ok = F)
  checkmate::assertNumeric(criticalValue, len = 1, null.ok = F)
  checkmate::assertNumeric(precision, len = 1, null.ok = F)
  checkmate::assertLogical(mean, len = 1, null.ok = F)
  if (is.null(outliers))
    joutliers<-.jnull("[Ljava/lang/String;")
  else
    joutliers=.jarray(outliers, "java.lang.String")
  jrslt<-.jcall("demetra/highfreq/r/FractionalAirlineProcessor", "Ljdplus/highfreq/FractionalAirlineEstimation;", "estimate", 
                as.numeric(y), 
                rjd3toolkit:::matrix_r2jd(x), mean, .jarray(periods), as.integer(ndiff), ar, joutliers
                , criticalValue, precision, approximateHessian)
  model<-list(
    y=as.numeric(y),
    variables=rjd3toolkit:::proc_vector(jrslt, "variables"),
    X=rjd3toolkit:::proc_matrix(jrslt, "regressors"),
    b=rjd3toolkit:::proc_vector(jrslt, "b"),
    bcov=rjd3toolkit:::proc_matrix(jrslt, "bvar"),
    linearized=rjd3toolkit:::proc_vector(jrslt, "lin")
  )
  estimation<-list(
    parameters=rjd3toolkit:::proc_vector(jrslt, "parameters"),
    score=rjd3toolkit:::proc_vector(jrslt, "score"),
    covariance=rjd3toolkit:::proc_matrix(jrslt, "pcov")
  )
  likelihood<-rjd3toolkit:::proc_likelihood(jrslt, "likelihood.")
  
  return(structure(list(
    model=model,
    estimation=estimation,
    likelihood=likelihood),
    class="JDFractionalAirlineEstimation"))
  
}

#' Title
#'
#' @param y 
#' @param periods 
#' @param ndiff 
#' @param stde 
#' @param nbcasts 
#' @param nfcasts 
#'
#' @return
#' @export
#'
#' @examples
multiAirlineDecomposition.raw<-function(y, periods, ndiff=2, ar=F, stde=F, nbcasts=0, nfcasts=0){
  checkmate::assertNumeric(y, null.ok = F)
  
  jrslt<-.jcall("demetra/highfreq/r/FractionalAirlineProcessor", "Ljdplus/highfreq/FractionalAirlineDecomposition;", "decompose", as.numeric(y), 
                .jarray(periods), as.integer(ndiff), ar, stde, as.integer(nbcasts), as.integer(nfcasts))
  
  return (jrslt)
}

#' Title
#'
#' @param jdecomp 
#'
#' @return
#' @export
#'
#' @examples
multiAirlineDecomposition.ssf<-function(jdecomp){
  jssf<-.jcall("demetra/highfreq/r/FractionalAirlineProcessor", "Ljdplus/ssf/extractors/SsfUcarimaEstimation;", "ssfDetails", jdecomp)
  return (new(Class= "JD3_ProcResults", internal=jssf))
}

#' Title
#'
#' @param y 
#' @param period 
#' @param sn 
#' @param stde 
#' @param nbcasts 
#' @param nfcasts 
#'
#' @return
#' @export
#'
#' @examples
fractionalAirlineDecomposition.raw<-function(y, period, sn=F, stde=F, nbcasts=0, nfcasts=0){
  checkmate::assertNumeric(y, null.ok = F)
  checkmate::assertNumeric(period, len = 1, null.ok = F)
  checkmate::assertLogical(sn, len = 1, null.ok = F)
  jrslt<-.jcall("demetra/highfreq/r/FractionalAirlineProcessor", "Ljdplus/highfreq/FractionalAirlineDecomposition;", "decompose", as.numeric(y), 
                period, sn, stde, as.integer(nbcasts), as.integer(nfcasts))
  return (jrslt)
}

#' Title
#'
#' @param jdecomp 
#'
#' @return
#' @export
#'
#' @examples
fractionalAirlineDecomposition.ssf<-function(jdecomp){
  jssf<-.jcall("demetra/highfreq/r/FractionalAirlineProcessor", "Ljdplus/ssf/extractors/SsfUcarimaEstimation;", "ssfDetails", jdecomp)
  return (new(Class= "JD3_ProcResults", internal=jssf))
}


#' Title
#'
#' @param jrslt 
#' @param stde 
#'
#' @return
#' @export
#'
#' @examples
jd2r_multiAirlineDecomposition<-function(jrslt, stde=F){
  
  #ucarima model
  ncmps<-rjd3toolkit:::proc_int(jrslt, "ucarima.size")
  model<-arima_extract(jrslt, "ucarima.model")
  cmps<-lapply(1:ncmps, function(cmp){return (ucm_extract(jrslt, cmp))})
  ucarima<-rjd3modelling::ucarima.model(model, cmps)
  
  yc<-rjd3toolkit:::proc_vector(jrslt, "y")
  estimation<-list(
    parameters=rjd3toolkit:::proc_vector(jrslt, "parameters"),
    score=rjd3toolkit:::proc_vector(jrslt, "score"),
    covariance=rjd3toolkit:::proc_matrix(jrslt, "pcov")
  )
  likelihood<-rjd3toolkit:::proc_likelihood(jrslt, "likelihood.")
  ncmps<-rjd3toolkit:::proc_int(jrslt, "ncmps")
  if (stde){
    decomposition<-lapply((1:ncmps), function(j){return (cbind(rjd3toolkit:::proc_vector(jrslt, paste0("cmp(",j, ")")),
                                                               rjd3toolkit:::proc_vector(jrslt, paste0("cmp_stde(",j, ")"))  ))})
  }else{
    decomposition<-lapply((1:ncmps), function(j){return (rjd3toolkit:::proc_vector(jrslt, paste0("cmp(",j, ")")))})
  }
  
  return(structure(list(
    ucarima=ucarima,
    decomposition=decomposition,
    estimation=estimation,
    likelihood=likelihood),
    class="JDFractionalAirlineDecomposition"))
}


#' Title
#'
#' @param jrslt 
#' @param sn 
#' @param stde 
#'
#' @return
#' @export
#'
#' @examples
jd2r_fractionalAirlineDecomposition<-function(jrslt, sn=F, stde=F){
  #ucarima model
  ncmps<-rjd3toolkit:::proc_int(jrslt, "ucarima.size")
  model<-arima_extract(jrslt, "ucarima.model")
  cmps<-lapply(1:ncmps, function(cmp){return (ucm_extract(jrslt, cmp))})
  ucarima<-rjd3modelling::ucarima.model(model, cmps)
  
  yc<-rjd3toolkit:::proc_vector(jrslt, "y")
  sa<-rjd3toolkit:::proc_vector(jrslt, "sa")
  s<-rjd3toolkit:::proc_vector(jrslt, "s")
  if (sn){
    if (stde){
      decomposition<-list(
        y=yc,
        sa=sa,
        s=s,
        s.stde=rjd3toolkit:::proc_vector(jrslt, "s_stde")
      )
    }else{
      decomposition<-list(
        y=yc,
        sa=sa,
        s=s
      )
    }
  }else{
    t<-rjd3toolkit:::proc_vector(jrslt, "t")
    i<-rjd3toolkit:::proc_vector(jrslt, "i")
    if (stde){
      decomposition<-list(
        y=yc,
        t=t,
        sa=sa,
        s=s,
        i=i,
        t.stde=rjd3toolkit:::proc_vector(jrslt, "t_stde"),
        s.stde=rjd3toolkit:::proc_vector(jrslt, "s_stde"),
        i.stde=rjd3toolkit:::proc_vector(jrslt, "i_stde")
      )
    }else{
      decomposition<-list(
        y=yc,
        t=t,
        sa=sa,
        s=s,
        i=i
      )
      
    }
  }
  estimation<-list(
    parameters=rjd3toolkit:::proc_vector(jrslt, "parameters"),
    score=rjd3toolkit:::proc_vector(jrslt, "score"),
    covariance=rjd3toolkit:::proc_matrix(jrslt, "pcov")
  )
  likelihood<-rjd3toolkit:::proc_likelihood(jrslt, "likelihood.")
  
  return(structure(list(
    ucarima=ucarima,
    decomposition=decomposition,
    estimation=estimation,
    likelihood=likelihood),
    class="JDFractionalAirlineDecomposition"))
}


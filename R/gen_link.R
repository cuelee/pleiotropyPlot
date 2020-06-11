#' gen_link
#'
#' This function help input data parsing
#'
#' @param factors example
#' @return A matrix of the infile
#' @export
gen_link = function(m, h2){
  f = rownames(m)

  fs = (apply(abs(m), MARGIN = 1, sum)-1)
  h2_mod = h2/max(h2)*0.8
  low = -1*h2_mod
  high = h2_mod

  bin = start = end = sign = NULL;
  for (bs in f){
    s = low[bs]
    for (ps in f){
      if(bs==ps) {next}
      else{
        e = high[bs]
        si = m[bs,ps]
        bin = c(bin,paste(bs,ps,sep='-'))
        start = c(start, s)
        end = c(end, e)
        sign = c(sign,si)
      }
    }
  }
  p = cbind(start,end,sign)
  row.names(p) =bin
  remove(bin,start,end,e,s,bs,ps,sign,si)

  lm = done = NULL;
  for (ffac in factors){
    for (tfac in factors[factors %notin% done]){
      if (ffac == tfac ) {next}
      else{
        lm = rbind(lm, c(ffac, p[paste(ffac,tfac,sep='-'),c('start','end','sign')], tfac, p[paste(tfac,ffac,sep='-'),c('start','end','sign')]))
      }
    }
    done = c(done,ffac)
  }
  colnames(lm) = c('f_factor','f_rstart','f_rend','f_sign','t_factor','t_rstart','t_rend','t_sign')
  ldf = data.frame(lm)

  return(ldf)
}

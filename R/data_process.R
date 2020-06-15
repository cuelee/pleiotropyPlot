#' gen_link
#'
#' This function help input data parsing
#'
#' @param dat$traits example
#' @return A matrix of the infile
gen_link = function(m, h2){
  `%notin%` <- Negate(`%in%`)
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
  for (ffac in f){
    for (tfac in f[f %notin% done]){
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


corrMatOrder <- function(
  corr,
  order = c("AOE", "FPC", "hclust", "alphabet"),
  hclust.method = c("complete", "ward", "ward.D", "ward.D2", "single",
                    "average", "mcquitty", "median", "centroid") )
{
  order <- match.arg(order)
  hclust.method <- match.arg(hclust.method)

  switch(order,
         AOE = reorder_using_aoe(corr),
         FPC = reorder_using_fpc(corr),
         hclust = reorder_using_hclust(corr, hclust.method),
         alphabet = sort(rownames(corr))
  )
}

#' Reorder the variables using the angular order of the eigenvectors.
reorder_using_aoe <- function(corr) {
  x.eigen <- eigen(corr)$vectors[, 1:2]
  e1 <- x.eigen[, 1]
  e2 <- x.eigen[, 2]
  alpha <- ifelse(e1 > 0, atan(e2 / e1), atan(e2 / e1) + pi)
  order(alpha) # returned vector
}

#' Reorder the variables using the first principal component.
reorder_using_fpc <- function(corr) {
  x.eigen <- eigen(corr)$vectors[, 1:2]
  e1 <- x.eigen[, 1]
  order(e1) # returned vector
}

#' Reorder the variables using hclust (Hierarchical Clustering).
reorder_using_hclust <- function(corr, hclust.method) {
  hc <- hclust(as.dist(1 - corr), method = hclust.method)
  order.dendrogram(as.dendrogram(hc)) # returned vector
}

maxp = function(x,thres = 30){
  return(max(10^(-thres),x))
}


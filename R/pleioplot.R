rainbow_palete_hex = c('#F44336','#2196F3','#8BC34A','#FF5722','#9C27B0','#00BCD4','#FFEB3B','#9E9E9E','#3F51B5','#4CAF50','#FF9800','#E91E63','#03A9F4','#CDDC39','#795548','#673AB7','#009688','#FFC107','#607D8B')
default_alpha = format(as.hexmode(floor(0.7*255)),width = 2, upper.case=T)

#' @export
EHtoSH = function(hex_col){
  white = '#FFFFFF'
  out_col = function(f){
    format(as.hexmode(floor(f)),width=2,upper.case=T)
  }
  R3 = as.hexmode(substr(hex_col,2,3))
  G3 = as.hexmode(substr(hex_col,4,5))
  B3 = as.hexmode(substr(hex_col,6,7))
  A1 = as.hexmode(substr(hex_col,8,9))
  R2 = as.hexmode(substr(white,2,3))
  G2 = as.hexmode(substr(white,4,5))
  B2 = as.hexmode(substr(white,6,7))
  #if(A1 == '0'){return(toupper(paste('#',R3,G3,B3,sep='')))}
  if(A1 == as.hexmode('ff')){return(toupper(paste('#',R2,G2,B2,sep='')))}
  op = (as.integer(A1))/255
  R1 = (as.integer(R3)*(op) + as.integer(R2)*(1-op))
  G1 = (as.integer(G3)*(op) + as.integer(G2)*(1-op))
  B1 = (as.integer(B3)*(op) + as.integer(B2)*(1-op))
  return(toupper( paste('#',out_col(R1), out_col(G1),  out_col(B1),sep='')))
}

#' pleioplot
#'
#' This function can generate circo plot.
#'
#' @param factors example
#' @return A matrix of the infile
#' @export
pleioplot =  function(snp, traits, rg_matrix, sumstats, pleioin, pleiores, h2, snp_reference, eta_col = c('#E77E23','#309F86'), rg_col = c('#F51929','#2341F6'), link_hex = 3, size_scale = 0.6, pleioplot_palette = paste(rainbow_palete_hex, default_alpha, sep='') ){

  ## gen inputs
  pleiores$pleio_p <- sapply(pleiores$pleio_p,maxp)
  traits = colnames(rg_matrix)
  pleio_out = pleiores[snp,]

  bp = dat$ref[snp,'BP']
  chr = dat$ref[snp,'CHR']

  bp_s = bp-1000000
  bp_e = bp+1000000

  ref_subset = dat$ref[which(dat$ref$CHR == chr & dat$ref$BP > bp_s & dat$ref$BP < bp_e),]
  ref_window_snps = rownames(ref_subset)[order(ref_subset$BP)]
  manhattan_snps = ref_window_snps[ref_window_snps %in% rownames(sumstats)]
  manhattan_p = sumstats[manhattan_snps, (1:length(traits))*2]

  manhattanplot_x = (ref_subset[manhattan_snps,'BP'] - bp) / 1000000 * 0.8
  manhattanplot_y = y_snp = NULL
  for(tr in traits){
    y_snp = c(y_snp, list( -log10(sumstats[snp,paste(tr,'_p',sep='')]) ))
    manhattanplot_y = c(manhattanplot_y, list(-log10(manhattan_p[manhattan_snps, paste(tr,'_p',sep='')])))
  }
  names(y_snp) = names(manhattanplot_y) = traits

  ldf = gen_link(rg_matrix, h2)

  sumstats_beta = sumstats[snp,paste(traits,'_beta',sep='')]
  names(sumstats_beta) = traits
  eta = pleioin[snp,paste(traits,'_beta',sep='')]
  names(eta) = traits
  eta_se = pleioin[snp,paste(traits,'_se',sep='')]
  names(eta_se)=traits

  lci = eta - eta_se * 1.96
  uci = eta + eta_se * 1.96
  eta_ci = list(lci, uci)

  col_mat = matrix(pleioplot_palette[1:length(traits)], nrow = 1)
  colnames(col_mat) = traits

  p=rep(0,length(traits))
  for (i in 1:length(traits)){
    p[i] = 10^(-y_snp[[traits[i]]])
  }
  names(p) = traits

  circlize::circos.clear()

  circlize::circos.par(start.degree = 0, canvas.xlim = c(-1/size_scale,1/size_scale), gap.degree = 2, cell.padding = c(0,0,0,0), canvas.ylim = c(-1.5,1.5), points.overflow.warning=F)
  circlize::circos.initialize(factors = traits, xlim = c(-1, 1))

  ##  Track 1
  cex_text = 0.6
  circlize::circos.track(factors = traits, x = rep(-0.3,length(traits)), y = rep(0, length(traits)), ylim = c(-0.01,0.01),
               bg.col = NA,
               bg.border = NA,
               track.height = 0.05,
               panel.fun = function(x, y) {
                 circlize::circos.text(x = x, y = circlize::CELL_META$cell.ylim[2] + circlize::uy(2, "mm"),
                             labels = sapply(circlize::CELL_META$sector.index, function(x) paste('',x,sep='')), facing='clockwise', niceFacing =T, adj = c(0,0), cex = cex_text, font = 4 )
                 circlize::circos.text(x = x+0.4, y = circlize::CELL_META$cell.ylim[2] + circlize::uy(2, "mm"),
                             labels = paste('P =', formatC(p[circlize::CELL_META$sector.index],format='E',digit=2),sep='') , facing='clockwise', niceFacing =T,  adj = c(0,0), cex = cex_text, font = 3)

                 circlize::circos.text(x = x+0.8, y = circlize::CELL_META$cell.ylim[2] + circlize::uy(2, "mm"),
                             labels = paste('BETA: ',formatC(as.numeric(sumstats_beta[circlize::CELL_META$sector.index]),format='E',digit=1,drop0trailing=T), sep = ''), facing='clockwise', niceFacing =T,  adj = c(0,0), cex = cex_text, font = 3)

                 circlize::circos.rect(xleft = -1, xright = 1, ybottom = -0.005, ytop = 0.005, border = NA, col= col_mat[1,circlize::CELL_META$sector.index])
               })

  ##  Track 2
  ymin = 0
  ymax = 0
  for(trait in traits){
    ymax = max(c(ymax, manhattanplot_y[[trait]]))
  }
  circlize::circos.track(ylim = c(ymin, ymax),
               x = rep(0,length(traits)),
               y = rep(0,length(traits)),
               bg.col = NA,
               bg.border = col_mat[, traits],
               bg.lty= 4,
               bg.lwd= 0.5,
               track.height = 0.2,
               panel.fun = function(x, y)
               {
                 circlize::circos.segments(x0=-0.9, y0 = 4, x1=0.9, y1=4, lty = 4, lwd = 0.1, col = '#FF000088')
               }
  )

  for (trait in traits){
    circlize::circos.points(x = manhattanplot_x,
                  y = manhattanplot_y[[trait]],
                  sector.index = trait,
                  track.index = 2,
                  pch = '.',
                  cex = 0.4,
                  #lwd = 0.3,
                  col = col_mat[, trait]
    )
    circlize::circos.points(x = 0,
                  y = y_snp[[trait]],
                  sector.index = trait,
                  track.index = 2,
                  pch = 20,
                  cex = 0.4,
                  col = '#000000AA'
    )
  }

  ##  Track 3
  yl = 1
  yp = 0.4
  ci_lwd=0.5
  circlize::circos.par("track.height" = 0.08)
  xl = eta_ci[[1]]/max(abs(eta))
  xu = eta_ci[[2]]/max(abs(eta))
  circlize::circos.track(ylim = c(-yl, yl),
               x = 0.8 * eta/(max(abs(eta))),
               y=rep(0,length(eta)),
               bg.col = NA,
               bg.border = col_mat[, traits],
               bg.lty= 4,
               bg.lwd= 0.5,
               track.height = 0.08,
               panel.fun = function(x, y)
               {
                 circlize::circos.segments(x0=0, y0 = -yl*0.8, x1=0, y1=yl*0.8, lty = 4, lwd = 0.5, col = '#000000AA')
                 if (x>0){
                   circlize::circos.rect(xleft = min(x,0), xright = max(x,0), ybottom = -yl*yp, ytop = yl*yp, border = NA, col= eta_col[1])
                 }
                 else if (x<0){
                   circlize::circos.rect(xleft = min(x,0), xright = max(x,0), ybottom = -yl*yp, ytop = yl*yp, border = NA, col= eta_col[2])
                 }

                 circlize::circos.segments(x0=max(-1, xl[1,circlize::CELL_META$sector.index]), y0=0, x1 = min(1, xu[1,circlize::CELL_META$sector.index]), y1=0,lwd = ci_lwd, col = '#000000AA')
                 circlize::circos.segments(x0=max(-1, xl[1,circlize::CELL_META$sector.index]), y0 = - yl*yp*0.4, x1=max(-1, xl[1,circlize::CELL_META$sector.index]), y1 = yl*yp*0.6,lwd = ci_lwd, col = '#000000AA')
                 circlize::circos.segments(x0= min(1,xu[1,circlize::CELL_META$sector.index]), y0 = - yl*yp*0.4, x1= min(1,xu[1,circlize::CELL_META$sector.index]), y1 = yl*yp*0.6,lwd = ci_lwd, col = '#000000AA')
               })

  ##  Track links
  hr = 0.5
  for(i in 1:nrow(ldf)){
    sign = as.numeric(as.matrix(ldf[i,'f_sign']))
    sector.index1=as.character(ldf[i,'f_factor'])
    point1=as.numeric(as.matrix(ldf[i,c('f_rstart','f_rend')]))
    sector.index2=as.character(ldf[i,'t_factor'])
    point2=as.numeric(as.matrix(ldf[i,c('t_rstart','t_rend')]))
    hex = format(as.hexmode(max(floor((abs(min(sign,1)))^(link_hex)*255-1),0)),width=2, upper.case=T)
    if(sign>0){
      circlize::circos.link(sector.index1, point1, sector.index2, point2, col=paste(rg_col[1], hex,sep=''), h.ratio=hr)
    }
    if(sign<0){
      circlize::circos.link(sector.index1, point1, sector.index2, point2, col=paste(rg_col[2], hex,sep=''), h.ratio=hr)
    }
  }
  circlize::circos.clear()

  ## Legend
  x_legend = rev(seq(-1,1,0.1))
  y_legend = c((paste(rg_col[1],format(as.hexmode(sapply(floor(abs(x_legend[x_legend>0])^link_hex*255-1), function(x) max(x,0))),width=2,upper.case=T),sep='')),'#FFFFFFFF', rev(paste(rg_col[2],format(as.hexmode(sapply(floor(abs(x_legend[x_legend>0])^link_hex*255-1),function(x) max(x,0))),width=2, upper.case=T),sep='')))

  col_fun_link = circlize::colorRamp2(x_legend, sapply(y_legend, EHtoSH))

  lgd_links = ComplexHeatmap::Legend(at = rev(seq(-1,1,0.5)), labels = rev(seq(-1,1,0.5)), col_fun = col_fun_link, title_position = "topleft", title = bquote(r[g]),labels_gp = grid::gpar(fontsize = 5), direction = "vertical", grid_width = grid::unit(0.3,'cm'))

  #draw(lgd_links, x = unit(0.03, "npc"), y = unit(0.2, "npc"), just = c("left", "top"))

  col_fun_beta = function(xs){
    res = rep(0,length(xs))
    for (i in 1:length(xs)){
      if (xs[i] >= 0) {
        res[i] = eta_col[1]
      }
      else if(xs[i]<0){
        res[i] = eta_col[2]}
    }
    return(res)
  }

  x_legend = eta/max(abs(eta))*0.8

  if(min(eta)>0){
    at_eta = c(formatC(as.numeric(max(eta)),format='E',digit=2) , 0)
  } else if(max(eta)<0){
    at_eta = c(formatC(as.numeric(min(eta)),format='E',digit=2), 0)
  } else {
    at_eta = c(formatC(as.numeric(min(eta)),format='E',digit=2), 0,formatC(as.numeric(max(eta)),format='E',digit=2))
  }

  lgd_eta = ComplexHeatmap::Legend(at = round(as.numeric(at_eta),3), labels = at_eta, col_fun = col_fun_beta, title = bquote(eta), labels_gp = grid::gpar(fontsize = 5), grid_width = grid::unit(0.3,'cm'))

  lgd_list_horizontal = ComplexHeatmap::packLegend(lgd_links, lgd_eta, direction = "horizontal")

  ComplexHeatmap::draw(lgd_list_horizontal, x = grid::unit(0.09, "npc"), y = grid::unit(0.99, "npc") , just = c("top"))
}

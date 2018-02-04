one.sample.ttest = function(m,s,mu,n){
    tstat = (m-mu)/s
    2*pt(abs(tstat), n-1, lower.tail=FALSE)
}


thermo_model_pred = function(s1,s2,x,y,B){

  X = (B+x)/(B+1)
  Y = (B+y)/(B+1)
  
  tmp = ((X-s1)*(Y-s1)*(X-s2)*(Y-s2)) / ((s1*(B+1)-B)*(s2*(B+1)-B)*(X-1)*(Y-1))
  b = -1*(X+Y+(B+1)*tmp)
  pr = 0.5 * (- b - sqrt(b^2-4*(B*tmp+X*Y)))
  
  # prediction is impossible if one of the two single mutants is out of the range
  # assign these predictions the max or min possible values (i.e. the asymptotes of the model)
  # if one single is above max and the other below min, the one that destroys the interaction has precedence
  max.ppi = min(c(X,Y))
  min.ppi = B/(B+1)
  pr[s1 > max.ppi | s2 > max.ppi] = max.ppi
  pr[s1 < min.ppi | s2 < min.ppi] = min.ppi
  
  pr
  
}



# function arguments for specific panels in the "pairs" plotting function
panel.hist40 <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks=40)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
panel.cor.pearson <- function(x, y, digits = 2, prefix = "", cex.cor,  ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, use="pairwise.complete.obs", method="pearson")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}
panel.smooth <- function(...){
  smoothScatter(..., nrpoints = 0, colramp = colorRampPalette(c("white", blues9[2:9])), add = TRUE)
}



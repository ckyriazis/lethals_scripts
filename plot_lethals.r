
### Script for plotting analytical predictions of MSB model
### and comparing against segregating recessive lethals
### for humans and Drosophila
# Author: Chris Kyriazis

# function for MSB predictions based on Nei 1968
calc_lethals <- function(u,N,coding_len,ns_ratio,frac){
  s=1
  q = u*sqrt(2*pi*N/s)
  L=coding_len*ns_ratio*frac*2

  num_lethals <- q*L
  
  return(num_lethals)
}


setwd("~/Documents/UCLA/lethals/analysis/plots/")
pdf("lethals_Ne.pdf", width = 10, height=5.5)

par(mfrow=c(1,2),  bty = "n")

width=4


# humans
lethal_frac <- seq(from = 0.001, to = 0.05, by = 0.001)
human_lethal_perc_N1 <- calc_lethals(u=1.5e-8,N=20000,coding_len= 30150000,ns_ratio = 0.7,frac = lethal_frac)
human_lethal_perc_N2 <- calc_lethals(u=1.5e-8,N=30000,coding_len= 30150000,ns_ratio = 0.7,frac = lethal_frac)
human_lethal_perc_N3 <- calc_lethals(u=1.5e-8,N=60000,coding_len= 30150000,ns_ratio = 0.7,frac = lethal_frac)

plot(lethal_frac*100,human_lethal_perc_N1, type = "l", xlab = "% lethal mutations", ylab = "# recessive lethals per diploid", lwd = width, main = "Humans", ylim=c(0,20), cex.lab=1.3, cex.axis=1.2)
lines(lethal_frac*100,human_lethal_perc_N2, lwd=width)
lines(lethal_frac*100,human_lethal_perc_N3, lwd=width)
x <- c(0,5,5,0)
y <- c(0.6,0.6,1.6,1.6)
polygon(x, y, col =rgb(1, 0, 0,0.4), border = F) 


# drosophila
lethal_frac <- seq(from = 0.001, to = 0.05, by = 0.001)
dmel_lethal_perc_N1 <- calc_lethals(u=3e-9,N=500000,coding_len= 22005320,ns_ratio = 0.74,frac = lethal_frac)
dmel_lethal_perc_N2 <- calc_lethals(u=3e-9,N=1000000,coding_len= 22005320,ns_ratio = 0.74,frac = lethal_frac)
dmel_lethal_perc_N3 <- calc_lethals(u=3e-9,N=5000000,coding_len= 22005320,ns_ratio = 0.74,frac = lethal_frac)

plot(lethal_frac*100,dmel_lethal_perc_N1, type = "l", xlab = "% lethal mutations", ylab = "# recessive lethals per diploid", lwd = width, main = "Drosophila", ylim=c(0,30), cex.lab=1.3, cex.axis=1.2)
lines(lethal_frac*100,dmel_lethal_perc_N2, lwd=width)
lines(lethal_frac*100,dmel_lethal_perc_N3, lwd=width)
x <- c(0,5,5,0)
y <- c(1,1,3,3)
polygon(x, y, col =rgb(1, 0, 0,0.4), border = F) 


dev.off()





# number of new lethals per diploid
#humans
2*1.5e-8*30000000*2.31/3.31*0.005
2*3e-9*22005320*0.74*0.01






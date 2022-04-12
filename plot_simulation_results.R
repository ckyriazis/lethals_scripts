
### script for plotting output of lethals SLiM simulation
### and add analytical prediction to plot
# Author: Chris Kyriazis


### read in data and plot
setwd("~/Documents/UCLA/lethals/analysis/data/")
data <- read.csv("lethals_sim_output.csv")


setwd("~/Documents/UCLA/lethals/analysis/plots//")
pdf("human_lethals_sim_plot.pdf",width = 8, height = 8)

xmin<- 0
xmax <- 16000

par(mfrow=c(2,1),mar = c(4,5,2,2), bty = "n")
plot(data$gen,data$popSize_AF, type="l", xlim=c(xmin,xmax),ylim=c(0,100000), lwd=4, xlab="Generation", ylab="N", col = "#D55E00", cex.lab=1.3)
lines(data$gen,data$popSize_EU, col="#56B4E9", lwd=4)
lines(data$gen,data$popSize_AS, col="#009E73", lwd=4)

legend(x="topleft",legend = c("Africa", "Europe","Asia"), col = c("#D55E00","#56B4E9","#009E73"), lty=1, lwd=4)

plot(data$gen,data$lethals_AF, type="l", xlim=c(xmin,xmax), ylim=c(0,2), lwd=2, xlab="Generation", ylab="# lethals per diploid", 
     col = "#D55E00", cex.lab=1.3)
lines(data$gen,data$lethals_EU, col="#56B4E9", lwd=2)
lines(data$gen,data$lethals_AS, col="#009E73", lwd=2)



### add analytical prediction to plot 
# function for MSB predictions based on Nei 1968
calc_lethals <- function(u,N,coding_len,ns_ratio,frac){
  s=1
  q = u*sqrt(2*pi*N/s)
  L=coding_len*ns_ratio*frac*2
  num_lethals <- q*L
  return(num_lethals)
}


lethals_AF <- c()
lethals_EU <- c()
lethals_AS <- c()


mu=1.5e-8
L=30158040
lethal_frac=0.005

for(i in 1:length(data$popSize_AF)){
  lethals_AF <- c(lethals_AF,calc_lethals(u=mu,N=data$popSize_AF[i],coding_len= L,ns_ratio = 0.7,frac = lethal_frac))
  lethals_EU <- c(lethals_EU,calc_lethals(u=mu,N=data$popSize_EU[i],coding_len= L,ns_ratio = 0.7,frac = lethal_frac))
  lethals_AS <- c(lethals_AS,calc_lethals(u=mu,N=data$popSize_AS[i],coding_len= L,ns_ratio = 0.7,frac = lethal_frac))
  
}
lines(data$gen,lethals_AF, lty=2, col = "#D55E00", lwd=3)
lines(data$gen,lethals_EU,lty=2, col = "#56B4E9", lwd=3)
lines(data$gen,lethals_AS,lty=2, col = "#009E73", lwd=3)


dev.off()











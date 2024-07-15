#!/usr/bin/env Rscript

ps = read.table("sizes.dat",header=F, sep=";")
pit = read.table("pores_in_time.dat", header=F, sep=";")
x = ps[-1,-1]
#x = cbind(ps,pit[,4])
x = as.matrix(x)
p = which(pit[,4]>10)
at=as.matrix(ps[1,-1])

cl = seq(1,length(p))
# time in ns (start, end, dt)
time = seq(0,999.5,0.5)

png("pore_in_time.png",width=1200,height=900)
for(i in cl){
    plot(time,x[p[i],], type="l", ylim=c(0,max(x)+0.2), col=i,cex.axis=2,cex.lab=2,lwd=3,xlab="time (ns)",ylab="Pore area nm^2")
    par(new=T)
}
legend("topleft",pch=15,col=cl, legend=pit[p,1])
dev.off()
mps = apply(x,2,sum)

pore_mps = NULL
for(i in seq(1,length(mps))){pore_mps=c(pore_mps,mps[i]/at[i]*100)}

png("percentage_pore.png",width=1200,height=900)
plot(time,pore_mps,type="l",cex.axis=2,cex.lab=2,lwd=3,xlab="time (ns)",ylab="Pore fraction %")
dev.off()

# files for xmgrace
mat = time
for (i in p){
    mat=cbind(mat,as.matrix(x[i,]))
}
write.table(mat, "pore_xvg.dat", sep = "\t", quote=F, col.names=F, row.names=F)


mat = time
mat=cbind(mat,as.matrix(pore_mps))
write.table(mat, "percentage_xvg.dat", sep = "\t", quote=F, col.names=F, row.names=F)

mat = time
mat=cbind(mat,as.matrix(mps))
write.table(mat, "pore_total_xvg.dat", sep = "\t", quote=F, col.names=F, row.names=F)


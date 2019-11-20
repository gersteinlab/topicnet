

argv=commandArgs(TRUE)
if (length(argv)==0 | length(argv)>10){
	print("usage: R -f lda.r --no-save --args mat num_kclass(4) output entropycutoff(0-1, recommend 0.1) var(2,3, how many trials/categories, var=1 if dist=bernoulli) dist(multinomial, bernoulli, rank ")
	q(save="no")
}

##argv=c("enhtf_tss.allbin.txt", "4", "enhtf_tss_bv1.rdata", "0", "1", "bernoulli")
#argv=c("enhtf_tss.allbin.txt", "4", "enhtf_tss_mv2.rdata", "0", "2", "multinomial")
#argv=c("enhtf_tss.allbin2.txt", "4", "enhtf_tss_m2v2.rdata", "0", "2", "multinomial")
#argv=c("enhtf_tss.allbin2.txt", "4", "enhtf_tss_m2v3.rdata", "0", "3", "multinomial")

input=argv[1]
K=as.integer(argv[2])
output=argv[3]
ecut=as.numeric(argv[4])
vvj=as.integer(argv[5])
dd=argv[6]

###will  vvj=1 dd="bernoulli"; vvj=2, dd="multinomial"; vvj=3, dd="multinomial"

print(paste("R -f lda.r --no-save --args ",input, K, output, ecut))

all.ds.bin=read.table(input, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE)
print(">>>Dimension of all.ds.bin:")
print(dim(all.ds.bin))

p.col.sumbin=colSums(all.ds.bin)/nrow(all.ds.bin)
ent.pcolsum=ifelse(p.col.sumbin==0 | p.col.sumbin==1, 0, -p.col.sumbin*log10(p.col.sumbin)-(1-p.col.sumbin)*log(1-p.col.sumbin))

keep.idx=which(ent.pcolsum>=ecut)
all.ds.bin=all.ds.bin[,keep.idx]
print(">>>Dimension of all.ds.bin after filtering:")
print(dim(all.ds.bin))

library("mixedMem")

#mm.list=list();
#ntimes=100
#for (xx in 1:ntimes){
 # mm.list[[xx]]=list()

 # cat("run ",xx,"\n")
  
  Total <- nrow(all.ds.bin)
                                        ## Number of variables
J <- ncol(all.ds.bin)
Rj <- rep(1, J)
Nijr <- array(1, dim = c(Total, J, max(Rj)))
#K <- 4

###0, or 1 two options
Vj <- rep(vvj, J)
alpha <- rep(.2, K)
dist <- rep(dd, J) ###

obs <- array(0, dim = c(Total, J, max(Rj), max(Nijr)))
obs[ , ,1,1] <- as.matrix(all.ds.bin)

ssdd=as.integer(runif(min=100,max=1000000,n=1))
set.seed(ssdd)
theta <- array(0, dim = c(J,K,max(Vj)))

if (dd=="bernoulli"){ ##bernoulli
for(j in 1:J)
{
  ##bernoulli dist
  theta[j,,] <- cbind(rbeta(K, 1,1), matrix(0, nrow = K, ncol = Vj[1]-1))

}
}else{
for(j in 1:J)
{
  ##multinomial dist
  theta[j, , ] <- gtools::rdirichlet(K, rep(.8, Vj[j]))
}

}

initial <- mixedMemModel(Total = Total, J = J, Rj = Rj, Nijr = Nijr, K = K, Vj = Vj, alpha = alpha, theta = theta, dist = dist, obs = obs)

out <- mmVarFit(initial)

#mm.list[[xx]][["seed"]]=ssdd
#mm.list[[xx]][["init"]]=initial
#mm.list[[xx]][["out"]]=out

save(list=ls(),file=paste(output,".k", K,".v",vvj,".d",dd, ".rdata",sep=""))


###add the output text and figures
###given the first half is gm, the left are k
##distribution of theta

ntf=nrow(all.ds.bin)
if (ntf %%2 ==0){
  ntf=ntf/2
  
  outf=output
  
  model=out
  lab.gm=rownames(all.ds.bin)[1:ntf]
  mem.est=model$phi/rowSums(model$phi)

  gm.est=mem.est[1:ntf,]
  k.est=mem.est[-c(1:ntf),]

  rownames(gm.est)=lab.gm
  rownames(k.est)=lab.gm

  ##class for each tf
  tf.cls=data.frame(gm=apply(gm.est, 1, function(x){ which.max(x);}), k=apply(k.est, 1, function(x){which.max(x);}) )


  
  ##plot

write.table(gm.est,file=paste(outf,"_k", K,"_v",vvj,"_d",dd,"_gm_est.txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T)
write.table(k.est,file=paste(outf,"_k", K,"_v",vvj,"_d",dd,"_k_est.txt",sep=""), sep="\t", quote=F, row.names=T, col.names=T)


  w=h=ceiling(sqrt(ntf))
  if (ntf%%w==0){
    h=h+1
  }


pdf(paste(outf,"_k", K,"_v",vvj,"_d",dd,"_diffclass4tf.pdf",sep=""), width=w*1.5, height=h*1.5)
  par(mfrow=c(w,h))
  #par(cex=8/K)
par(oma=c(3,5,3,1))
par(mar=rep(0.1,4))

indices=1:ntf

for(i in indices){
      plot(gm.est[i,], type="p", lwd=2, col="black", ylim=c(-0.1,1.1), xlim=c(0.25, (K+0.25)), yaxt="n", xaxt="n", pch=c(0:14))    
      text(0, 1, labels=lab.gm[i], cex=0.9, adj=c(0,1), pos=4)
      points(c(1:K), k.est[i,], col="red", pch=c(0:16))

}

plot(NA, type="n", xaxt="n", yaxt="n", ylim=c(-0.1,1.1), xlim=c(0.25, (K+.25)) )

#par(fig=c(0,1,0,1), oma=c(0,5, 0, 1), mar=rep(0,4), new=T)

legend("bottom", c("GM12878", "K562"), col=c("black","red"), cex=0.8,fill=c("black","red") )

dev.off()

print(paste(outf,".k", K,".v",vvj,".d",dd,"_diffclass4tf.pdf, has generated",sep=""))
  ########




###distribution  fo theta
theta.cls=out$theta[,,1]
if (dd=="multinomial"){
  theta.cls=out$theta[,,2]
}
rownames(theta.cls)=colnames(all.ds.bin)
w=ceiling(sqrt(K))
pdf(paste(outf,"_k", K,"_v",vvj,"_d",dd,"_mmsb_theta_dist.pdf",sep=""),width=w*3,height=w*3)
par(mfrow=c(w,w))
par(cex=0.7)
par(lend=2)
par(tcl= -0.15)   
par(mar=c(3,3,1,1)+0.1)
par(mgp=c(1.1, 0.15, 0))

  
for (i in 1:K){
 plot(density(theta.cls[,i]),xlab="Theta",ylab="Dist", main=paste("group",i))
}
dev.off()

print(paste(outf,"_k", K,"_v",vvj,"_d",dd,"_mmsb_theta_dist.pdf, generated",sep=""))
  

  
d_h=rowSums((sqrt(gm.est)-sqrt(k.est))**2)*(1/sqrt(2))

d_h.sort=sort(d_h, decreasing=T)

write.table(d_h.sort, file=paste(outf,"_k", K,"_v",vvj,"_d",dd,"_mm_rankdf.txt",sep=""), sep="\t",quote=F, row.names=T, col.names=T)

pdf(paste(outf,"_k", K,"_v",vvj,"_d",dd,"_mm_rank_diff.pdf",sep=""),width=12,height=6)
                               
xbar=barplot(d_h.sort,names.arg=rep("",length(d_h.sort)), col="lightblue", cex.axis=1,main="",ylim=c(-0.1,1.6),ylab="Hellinger Distance (GM versus K)",cex.lab=1)

text(xbar+0.4, d_h.sort+0.1, names(d_h.sort), srt=60,cex=0.7)

dev.off()  
print(paste(outf,"_k", K,"_v",vvj,"_d",dd,"_mm_rank_diff.pdf, generated",sep=""))

  
###difference of TF


}




#vizMem(out)

#enhtf_tss.allbin.txt", "4", "enhtf_tss_bv1.rdata", "0", "1", "bernoulli
#for n in `seq 1 100`; do for y in 4 8 16; do echo -e "#BSUB -W 96:00\n#BSUB -M 40000000\n#BSUB -q gerstein\n#BSUB -o raw_${n}_${y}.out\nmodule load Apps/R/3.1.1-generic\ncd `pwd` \n" >> final_${n}_${y}.sh;  for x in `find ~/scratch/mixmb/data/ -name "*.txt"` ;  do fn=`basename $x`;  echo -e "R -f ~/codes/lda.r --args $x $y ${fn}_${n}_$y_bv1 0 1 bernoulli \nR -f ~/codes/lda.r --args $x $y ${fn}_${n}_$y_mv2 0 2 multinomial  " >> final_${n}_${y}.sh; done; done; done


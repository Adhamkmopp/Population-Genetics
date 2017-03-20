
palette(c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac"))
palette(c("red","black","green","orange","blue","gray"))

##########################################################################
####################33
manPlot<-function(p,chr,mat,maxLogPval=10,...){

    keep<-!is.na(p)
    p<-p[keep]
    chr<-chr[keep]
    mat<-mat[keep,]
    p[-log(p,base=maxLogPval)>10]<-10^(-maxLogPval)

    bon <- -log(0.05/length(p),base=10)
#    cat("Bonferroni",bon,"\n")
    ylim<-c(0,max(c(-log(p,base=10),na.rm=T),bon)*1.05)
    plot(-log(p,base=10),col=ifelse(chr%%2==0,"#2166ac","#d1e5f0"),ylab=expression(-log[10](p-value)),xlab="Chromosomes",axes=F,pch=20,ylim=ylim, ...)
    axis(2)
    t<--log(0.05/nrow(p),base=10)
    abline(h=bon,lty=2,col=1,lwd=2)
    legend("topright","Bonferroni correction",lty=2,bty="n",lwd=2)
   
    ord<-order(p)

    cat("-------------- TOP 10 -----------------------\n")
    print(mat[head(ord,10),])
    
}
qqp<-
function(x,ci=TRUE,add=FALSE,ylab="Observed log10(p-value)",xlab="Expected log10(p-value)",maxLogP,...){
  x<-x[!is.na(x)]
  if(!missing(maxLogP))
    x[x<10^-maxLogP]<-10^-maxLogP
  N<-length(x)
 x<-sort(x)

#lambda<-round(median(x)/qchisq(0.5,1),2)
  e<- -log((1:N-0.5)/N,10)
  if(add)
    points(e,-log(x,10),...)
  else{
    plot(e,-log(x,10),ylab=ylab,xlab=xlab,...)
    abline(0,1,col=2,lwd=2)
  }
#legend("topleft",paste("lambda=",lambda))

  if(ci){
    c95<-qbeta(0.95,1:N,N-(1:N)+1)
    c05<-qbeta(0.05,1:N,N-(1:N)+1) 
    lines(e,-log(c95,10))
    lines(e,-log(c05,10))
  }
}


hmmtest<-function(name,name2) {
x<-scan(name2)
    pdfname<-paste(name,".png",sep="")
    cat("Figure name: ",pdfname)
bitmap(pdfname)
plot(scan(name),ylab="Statistic",xlab="position (SNP number)")
abline(h=quantile(x[-1],0.95),lty=2)
legend("topleft","significance threshold",lty=2)
dev.off()
}

post<-function(name,name2){
  cat(name2)
post<-read.table(name,na="-1")

  mat<-matrix(NA,ncol=20,nrow=20)
  print(dim(post))
mat[lower.tri(mat)]<-post[,as.integer(name2)]
pdfname<-paste(name,".png",sep="")
cat("Figure name: ",pdfname)
bitmap(pdfname)
print(lattice::levelplot(1-mat,ylab="ind", xlab="ind",main="IBD probabilities at a locus"))
dev.off()
}


acas<-function(name){
    pdfname<-paste(name,".png",sep="")
    cat("Figure name: ",pdfname)
bitmap(pdfname)
plot(scan(name),ylab="Fraction of cases belong to pop 2",xlab="SNP along a chromosomal region")
dev.off()
bitmap(sub("acas","test1",pdfname))
plot(scan(sub("acas","test1",name))[-c(1,201)],ylab="test statistic",xlab="SNP along a chromosomal region")
    abline(h=1,lty=2)
    legend("topleft","significance threshold",)
dev.off()
  }

pair<-function(name='plink.genome'){
#plot pairwise relatedness
pdfname<-paste(name,"pairwise_relatedness.pdf",sep="")
d<-read.table(name,hea=T,colC=c("character","character","character","character","character",rep("numeric",14-5)))

pdf(pdfname)
pair.pic(d,pdfname)
dev.off()
jpgname<-paste(name,"pairwise_relatedness",sep="")

pair.pic.jpeg(d,jpgname)
    cat("Pdf figure name: ",pdfname,"\n")
    cat("Jpg figure names: ",jpgname,"*.jpg\n")
}
pair.pic<-function(d,name){
  plot(d$Z1,d$Z2,ylab="k2",xlab="k1",main=" pairwise relatedness",col=6,lwd=2,ylim=c(0,1),xlim=c(0,1))

text(0,1,"MZ",font=2)
text(1,0,"PO",font=2)
text(0.5,0.25,"FS",font=2)
text(0.5,0,"HS",font=2)
text(0.25,0,"C1",font=2)
text(1/16,0,"C2",font=2)
text(0,0,"U",font=2)

plot(d$PI_HAT,ylab="relatedness (r)",xlab="pairs of individuals",main="relatedness",col="blue")

plot(0:1,0:1,col="transparent",ylab="",xlab="",axes=F)
ord<-order(d$PI_HAT,decreasing=T)

for(tal in 1:10){
g<-paste(d$FID1[ord[tal]]," vs ",d$FID2[ord[tal]])
k<-paste("k0=",d$Z0[ord[tal]]," k1=",d$Z1[ord[tal]]," k2=",d$Z2[ord[tal]]," r=",d$PI_HAT[ord[tal]]," p=",d$P[ord[tal]])
text(0.5,tal/10,g,cex=0.7,font=2)
text(0.5,tal/10+0.05,k,cex=0.7,font=2)
}

}
pair.pic.jpeg<-function(d,name){
  ##first
  png(paste(name,"_1.png",sep=""))
  plot(d$Z1,d$Z2,ylab="k2",xlab="k1",main=" pairwise relatedness",col=6,lwd=2,ylim=c(0,1),xlim=c(0,1))
  text(0,1,"MZ",font=2)
  text(1,0,"PO",font=2)
  text(0.5,0.25,"FS",font=2)
  text(0.5,0,"HS",font=2)
  text(0.25,0,"C1",font=2)
  text(1/16,0,"C2",font=2)
  text(0,0,"U",font=2)
  dev.off()
  ##second
  png(paste(name,"_2.png",sep=""))
  plot(d$PI_HAT,ylab="relatedness (r)",xlab="pairs of individuals",main="relatedness",col="blue")
  dev.off()
  ##third
  png(paste(name,"_3.png",sep=""),,height=540)
  plot(0:1,0:1,col="transparent",ylab="",xlab="",axes=F)
  ord<-order(d$PI_HAT,decreasing=T)

  for(tal in 1:10){
    g<-paste(d$FID1[ord[tal]]," vs ",d$FID2[ord[tal]])
    k<-paste("k0=",d$Z0[ord[tal]]," k1=",d$Z1[ord[tal]]," k2=",d$Z2[ord[tal]]," r=",d$PI_HAT[ord[tal]]," p=",d$P[ord[tal]])
    text(0.5,tal/10,g,cex=0.7,font=4)
    text(0.5,tal/10+0.05,k,cex=0.7,font=4)
  }
  dev.off()
}

###################
het<-function(name='plink.het'){

pdfname<-paste(name,"inbreeding.pdf",sep="")
d<-read.table(name,hea=T,colC=c("character","character",rep("numeric",4)))
pdf(pdfname)
#postscript("inbreeding.eps")
col="darkblue"#c(rep("darkblue",94-11),rep("darkred",11))
plot(d$F,lwd=1,ylab="Inbreeding (F)",xlab="individuals",main="",col=col) 
legend(5,0.18,c("Danish","Greenlanders"),col=c("darkblue","darkred"),bty="n",pch=1)

plot(0:1,0:1,col="transparent",ylab="",xlab="",axes=F)
ord<-order(d$F,decreasing=T)
for(tal in 1:10){
g<-paste(d$FID[ord[tal]], "F=",d$F[ord[tal]])
text(0.5,tal/10,g,cex=0.7,font=2)
}
dev.off()
}


###############
miss<-function(name='plink.lmiss'){

d<-read.table(name,hea=T,colC=c("NULL","NULL","integer","NULL"))
pdfname<-paste(name,".pdf",sep="")
pdf(pdfname)
hist(d[,1],br=100,col=1,main="Missingness",ylab="Number of SNPs",xlab="Number of missing genotypes")
dev.off()
}
########################
hwe<-function(name='plink.hwe'){

pdfname<-paste(name,".pdf",sep="")
h<-read.table(name,hea=T,nr=1)
c<-ifelse(names(h)=="P","numeric","NULL")
d<-read.table(name,hea=T,colC=c)
pdf(pdfname)
hist(d[,1],br=1000,col=2,main="HWE p-value",ylab="number of SNPs",xlab="p-value")
hist(d[,1],br=1000,col=2,main="HWE p-value",ylab="number of SNPs",xlab="p-value",xlim=c(0,0.9),ylim=c(0,10000))

dev.off()
}

sexcheck<-function(name='plink.sexcheck'){

pdfname<-paste(name,".pdf",sep="")
h<-read.table(name,hea=T,nr=1)
c<-ifelse(names(h)=="F","numeric","NULL")
c<-ifelse(names(h)=="PEDSEX","character",c)
c<-ifelse(names(h)=="FID","character",c)
c<-ifelse(names(h)=="STATUS","character",c)
d<-read.table(name,hea=T,colC=c)
pdf(pdfname)
 boxplot(d$F~d$PEDSEX,col=2:3,ylab="Inbreeding coefficient (X CHR)",names=c("males","females"),main=name)


plot(0:1,0:1,col="transparent",ylab="",xlab="",axes=F)

temp<-d[d$STATUS!="OK",]

for(tal in 1:nrow(temp)){
g<-paste(temp$FID[tal],"  F=",round(temp$F[tal],2))
text(0.5,tal/nrow(temp),g,cex=0.7,font=2)
}

dev.off()
}

##########################
freq<-function(name='plink.frq'){

pdfname<-paste(name,".pdf",sep="")
d<-read.table(name,hea=T,colC=c("NULL","NULL","NULL","NULL","numeric","NULL"))

pdf(pdfname)
hist(d[,1],br=1000,col=2,main="Frequency spectrum",ylab="number of SNPs",xlab="MAF")
hist(d[,1],br=1000,col=2,main="Frequency spectrum",ylab="number of SNPs",xlab="MAF",xlim=c(0.01,0.5),ylim=c(0,10000))

dev.off()

}
#############################
# name<-"data/plink.mibs"
# cc<-"data/gwa.fam"

cluster<-function(name="plink.mdist",cc=NA){

    status<-NA
  if(!is.na(cc)){
     id<- read.table(paste0(name,".id"))
      status<-read.table(cc,as.is=T)
      status<-status[status[,1]%in%id[,1],6]
      status[status==-9]<-3
  }
  
    
  m <- as.matrix(read.table(name))
  mds <- cmdscale(as.dist(1-m))
  pdfname<-paste(name,"cluster.pdf",sep="")
 
  pdf(pdfname)
  if(is.na(status[1]))
      plot(mds,col=2,lwd=2,ylab="PC2",xlab="PC1",main="Multidimensional scaling")
  else{
      plot(mds,lwd=2,ylab="Dimension 1",xlab="Dimension 2",main="Multidimensional scaling",col=c(1,2,4)[status])
      legend("topright",c("Controls","Cases","Unknown"),col=c(1,2,4),pch=1)
  }
  dev.off()
     cat("Figure name: ",pdfname,"\n")
}

#############################
# name<-"plink.mds"
# cc<-"data/gwa.fam"

mdsPlot<-function(name="plink.mds",cc=NA){

    pdfname<-paste(name,".pdf",sep="")
    pdf(pdfname)
    x = read.table(name,as.is=T,header=T)
    if(!is.na(cc)){
        y = read.table(cc,as.is=T)
        row.names(y)=y[,1]
        cols = y[x[,1],6]
    }else{
        cols = rep(4,nrow(x))
    }
    cols[is.na(cols)]=4
    cols[cols==-9]=4
    plot(x[,"C1"],x[,"C2"],col=cols,main="MDS",xlab="Dimension 1",ylab="Dimension 2")
    legend("topleft",c("Cases","Controls","Unknown disease status"),col=c(2,1,4),pch=21)
    dev.off()
    cat("Figure name: ",pdfname,"\n")
}



##########################################
qassoc.perm<-function(name="plink.qassoc.perm",...){
    figName<-paste(name,".png",sep="")
    bitmap(figName)
    p<-read.table(name,colC=c("integer",rep("NULL",2),"numeric","NULL"),head=T)
    plot(-log(p[,2],base=10),col=p[,1]%%2+3,...)
    t<--log(0.05/nrow(p),base=10)
    abline(h=t,lty=2,col=1)
    legend(0,t-0.5,"bonferroni correction",lty=2,bty="n")
    dev.off()
    cat("Figure:",figName,"\n")
}
##########################################
qassoc<-function(name="plink.qassoc",...){
    figName<-paste(name,".png",sep="")
    bitmap(figName)
    p<-read.table(name,colC=c("integer",rep("NULL",7),"numeric"),head=T)
    plot(-log(p[,2],base=10),col=p[,1]%%2+3,...)
    t<--log(0.05/nrow(p),base=10)
    abline(h=t,lty=2,col="black",lwd=2)
    legend(0,t-0.5,"Bonferroni correction",lty=2,bty="n",lwd=2)
    dev.off()
 cat("Manhattan plot can now be found in:",figName,"\n")

}
##########################################
logistic<-function(name="plink.assoc.logistic",ylim=NULL,...){
    figName<-paste(name,".png",sep="")
bitmap(figName)
p<-read.table(name,colC=c("integer","NULL","integer","NULL","character",rep("NULL",3),"numeric"),head=T)
p<-p[p$TEST=="ADD",]
if(is.null(ylim))
 ylim<-c(0,max(-log(p$P,base=10),7,na.rm=T))

plot(-log(p$P,base=10),col=p$CHR%%2+3,ylab=expression(-log[10](p-value)),ylim=ylim,pch=20,...)
t<--log(0.05/nrow(p),base=10)
abline(h=t,lty=2,col=1)
legend(0,t-0.5,"bonferroni correction",lty=2,bty="n")
dev.off()
cat("Figure:",figName,"\n")
} 
assoc<-function(name="plink.assoc",...){
    figName<-paste(name,".png",sep="")
bitmap(figName)
p<-read.table(name,colC=c("integer",rep("NULL",7),"numeric","NULL"),head=T)
plot(-log(p[,2],base=10),col=p[,1]%%2+3,ylab=expression(log[10](p-value)),xlab="Chromosomes",...)
t<--log(0.05/nrow(p),base=10)
abline(h=t,lty=2,col="black",lwd=2)
legend("topleft","Bonferroni correction",lty=2,bty="n",lwd=2)
dev.off()
cat("Manhattan plot can now be found in:",figName,"\n")

}
genomic.control<-function(name="plink.assoc",...){
figName<-    paste(name,".QQ.png",sep="")
bitmap(figName)

p<-read.table(name,colC=c("NULL",rep("NULL",6),"numeric","NULL","NULL"),head=T)
qqplot(rchisq(1000000,1),p$CHISQ,ylab="observed",xlab="Expected")
qqline(rchisq(1000000,1),col=1)
l<-median(p$CHISQ,na.rm=T)/0.456
text(2,max(p$CHISQ,na.rm=T),paste("lamda=",round(l,3)),bty="n") 
dev.off()
cat("QQ plot can now be found in:",figName,"\n")
}
genomic.controlQ<-function(name="plink.assoc",...){
figName<-    paste(name,".QQ.png",sep="")
bitmap(figName)
p<-read.table(name,head=T)
q<-qqp(p$P,pch=16,col="darkblue")
dev.off()
cat("Figure:",figName,"\n")
}
genomic.controlReg<-function(name="plink.assoc",...){
figName<-    paste(name,".QQ.png",sep="")
bitmap(figName)

p<-read.table(name,head=T)
p<-p[p$TEST=="ADD"|p$TEST=="REC"|p$TEST=="DOM" ,]
q<-qqp(p$P,pch=16,col="darkblue")
dev.off()
cat("Figure:",figName,"\n")
}
fisher<-function(name="plink.fisher",...){
bitmap(paste(name,".png",sep=""))
p<-read.table(name,colC=c("integer",rep("NULL",6),"numeric","NULL"),head=T)
plot(-log(p[,2],base=10),col=p[,1]%%2+3,...)
t<--log(0.05/nrow(p),base=10)
abline(h=t,lty=2,col=1)
legend(0,t-0.5,"bonferroni correction",lty=2,bty="n")
dev.off()
}


linear<-function(name="plink.assoc.linear",ylim=NULL,...){
    p<-read.table(name,colC=c("integer","character","integer","NULL","character",rep("NULL",3),"numeric"),head=T)
    p<-p[p$TEST=="ADD"|p$TEST=="DOM"|p$TEST=="REC",]

    bitmap(paste(name,".png",sep=""),res=200)
    manPlot(p$P,p$CHR,p)
    dev.off()

    
}



qassoc.mperm<-function(name="plink.qassoc.perm",ylim=c(0,1),...){
bitmap(paste(name,".png",sep=""))
p<-read.table(name,colC=c("integer",rep("NULL",2),"numeric","numeric"),head=T)

plot(p$EMP2,col=p$CHR%%2+3,ylim=ylim,ylab="corrected p-value",...)
t<-0.05
abline(h=t,lty=2,col=1)
legend(0,t-0.5,"permutation correction",lty=2,bty="n")
dev.off()
}
missing_fun<-function(name="plink.qassoc.perm",ylim=c(0,1),...){
  filename<-paste(name,".png",sep="")
 filename2<-paste(name,"2.png",sep="")
  cat("Figure name: ",filename,"\n")
  cat("Figure name: ",filename2,"\n")
bitmap(filename)
dat<-read.table(name,hea=T,colC=c("integer","NULL",rep("numeric",3)))
dat<-dat[!(dat[,2]==0&dat[,3]==0),]
hist(dat[,4],br=10,main="Test for differential missingness",xlab="p-value",col=1)
dev.off()
  bitmap(filename2)
  plot(-log(dat[,4],base=10),col=c(3,4)[dat[,1]%%2+1],ylab=expression(log[10](p-value)),...)
  print(-log10(0.05/sum(!is.na(dat[,4]))))
  abline(h=-log10(0.05/sum(!is.na(dat[,4]))),lty=2,col="black")
  legend("topleft","Bonferroni correction",lty=2,bty="n",lwd=2)
  dev.off()
}

##############################
tal<-1
l<-commandArgs(TRUE)
print(l)
t<-unlist(strsplit(l[tal],".",fixed=T))
if(length(t)==1)
    t<-unlist(strsplit(l[tal],"_",fixed=T))
type<-t[length(t)]
cat("type=",type,"\n")
{
if(type=="genome")
try(pair(l[tal]))
else if(type=="het")
try(het(l[tal]))
else if(type=="lmiss")
try(miss(l[tal]))
else if(type=="missing")
try(missing_fun(l[tal]))
else if(type=="hwe")
try(hwe(l[tal]))
else if(type=="sexcheck")
try(sexcheck(l[tal]))
else if(type=="frq")
try(freq(l[tal]))
else if(type=="mdist"|type=="mibs")
try(cluster(l[tal],l[tal+1]))
else if(type=="mds")
try(mdsPlot(l[tal],l[tal+1]))
else if(type=="perm")
try(qassoc.perm(l[tal]))
else if(type=="mperm")
try(qassoc.mperm(l[tal]))
else if(type=="qassoc"){
    try(qassoc(l[tal]))
    try(genomic.controlQ(l[tal]))
}
else if(type=="logistic"){
    try(logistic(l[tal]))
    try(genomic.controlReg(l[tal]))
}
else if(type=="linear"){
try(linear(l[tal]))
  try(genomic.controlReg(l[tal]))
}
else if(type=="assoc"){
try(assoc(l[tal]))
try(genomic.control(l[tal]))
}
else if(type=="fisher")
try(fisher(l[tal]))
else if(type=="acas")
try(acas(l[tal]))
else if(type=="txt")
try(hmmtest(l[tal],l[tal+1]))
else if(type=="post")
try(post(l[tal],l[tal+1]))


}
q("no")










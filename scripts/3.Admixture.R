library("snpStats")
library("Matrix")
library("survival")
library("RColorBrewer")
species <- read.plink('grants_filtered')
localities <- read.plink('localities')
snpk2=read.table("grants_filtered.2.Q")
snpk3=read.table("grants_filtered.3.Q")
snpk4=read.table("grants_filtered.4.Q")
snpk5=read.table("grants_filtered.5.Q")
df <-read.table(file='local', sep=" ")
locality <-df$V1
names=species$fam$pedigree
first<- snpk2[22,]
second <-snpk2[23,]
snpk2[22,] = snpk2[36,]
snpk2[23,] = snpk2[37,]
snpk2[36,] = first
snpk2[37,] = second
first<- snpk3[22,]
second <-snpk3[23,]
snpk3[22,] = snpk3[36,]
snpk3[23,] = snpk3[37,]
snpk3[36,] = first
snpk3[37,] = second
first<- snpk4[22,]
second <-snpk4[23,]
snpk4[22,] = snpk4[36,]
snpk4[23,] = snpk4[37,]
snpk4[36,] = first
snpk4[37,] = second
first<- snpk5[22,]
second <-snpk5[23,]
snpk5[22,] = snpk5[36,]
snpk5[23,] = snpk5[37,]
snpk5[36,] = first
snpk5[37,] = second
local_switch <- locality[22:23]
locality[22:23] <-locality[36:37]
locality[36:37] <-local_switch
graphics.off()
cols<-brewer.pal(n=6,name="Set1")
par(mfrow=c(4,1), xpd = T, mar = par()$mar + c(0,0,0,2))
par(fig=c(0, 1, 0, 0.25))
barplot(t(as.matrix(snpk2)), col= cols, border=NA,
        names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry, k=2")
legend(x=115, y=1,box.lwd=0.5,legend=c("granti","petersii"), fill = c(cols[2],cols[1]))
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , cex.axis=0.8, tck=FALSE, labels =locality)

par(fig=c(0, 1, 0.3, 0.55),new=TRUE)
barplot(t(as.matrix(snpk3)), col= cols, border=NA, names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry,k=3")
legend(x=115, y=1,legend=c("granti","petersii", 'notata'), fill = c(cols[2],cols[1], cols[3], cols[4]))
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , tck=FALSE, labels =locality)

par(fig=c(0, 1, 0.6, 0.85),new=TRUE)
barplot(t(as.matrix(snpk4)), col= c(cols[4],cols[2],cols[3],cols[1]),
        border=NA, names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry, k=4")
legend(x=115, y=1,legend=c("granti","petersii", "notata", "robertsii"), fill = c(cols[2],cols[1], cols[3],cols[4]))
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , tck=FALSE, labels =locality)
par(fig=c(0, 1, 0.3, 0.55),new=TRUE)

barplot(t(as.matrix(snpk5)), col= c(cols[3],cols[2],cols[5],cols[1], cols[4]),
        border=NA,  names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry, k=5")
legend(x=115, y=1,legend=c("granti","petersii", 'notata', 'thomson', 'robertsii'), fill = c(cols[2],cols[1],cols[3],cols[5],cols[4]))
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , cex.axis=0.8, tck=FALSE, labels =locality)

localities <- read.plink('localities')
snpk2=read.table("localities.2.Q")
snpk3=read.table("localities.3.Q")
snpk4=read.table("localities.4.Q")
snpk5=read.table("localities.5.Q")
snpk6=read.table("localities.6.Q")
snpk7=read.table("localities.7.Q")
snpk8=read.table("localities.8.Q")
names <-localities$fam$pedigree
species <-localities$fam$member
cols<-brewer.pal(n=8,name="Set1")
par(mfrow=c(4,1))
par(fig=c(0, 1, 0, 0.25))
barplot(t(as.matrix(snpk2)), col= cols, border=NA,
        names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry, K=2")
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , tck=FALSE, labels =species)

par(fig=c(0, 1, 0.3, 0.55),new=TRUE)
barplot(t(as.matrix(snpk3)), col= cols, border=NA,
        names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry, K=3")
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , tck=0, labels =species)

par(fig=c(0, 1, 0.6, 0.85),new=TRUE)
barplot(t(as.matrix(snpk4)), col= cols, border=NA,
        names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry, K=4")
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , tck=0, labels =species)

par(mfrow=c(4,1))
par(fig=c(0, 1, 0, 0.25))
barplot(t(as.matrix(snpk5)), col= cols, border=NA,
        names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry, K=5")
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , tck=0, labels =species)

par(fig=c(0, 1, 0.3, 0.55),new=TRUE)
barplot(t(as.matrix(snpk6)), col= cols, border=NA,
        names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry, K=6")
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , tck=0, labels =species)

par(fig=c(0, 1, 0.6, 0.85),new=TRUE)
barplot(t(as.matrix(snpk7)), col= cols, border=NA,
        names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry, K=7")
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , tck=0, labels =species)

par(mfrow=c(4,1))
par(fig=c(0, 1, 0.40, 0.65))
barplot(t(as.matrix(snpk8)), col= cols, border=NA,
        names.arg=(names), cex.names=0.8, las=2, ylab="Ancestry, K=8")
par(las=2)
axis(3, seq(0.5, 113.4,1.2) , tck=0, labels =species)

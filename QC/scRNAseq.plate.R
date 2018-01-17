# To run the script, put the count table and the index files in the same directory as the script is and in linux just type: Rscript scRNAseq.plate.R
# column names of the table should be in this form: "plateA_cell.index" in which the 'index' should be the number of the barcodes used per cell.

tab <- read.table("data.table.txt",header=T,sep="\t",row.names=1,fill=T)
tab <- na.omit(tab)
samples <- c("plateA","plateB","plateC")
pdf(file="quality_control_scRNAseq_384wplate.pdf",width=14,height=14)
par(mfrow=c(4,4))
for (j in samples){

umg <- tab[,grep(j,colnames(tab))]

T=colSums(umg)
ER=colSums(umg[grep("^ERCC",row.names(umg)),])
MT=colSums(umg[grep(":MT-",row.names(umg)),])
U=T-(ER+MT)
id <- read.table("./index.plate.txt")

comx <- NULL
for (i in 1:ncol(umg)) {comx <- c(comx,nrow(umg[umg[,i]>=1,]))}

ERcomx <- NULL
ERumg <- umg[grep("^ERCC",row.names(umg)),]
for (i in 1:ncol(ERumg)) {ERcomx <- c(ERcomx,nrow(ERumg[ERumg[,i]>=1,]))}

ss <- cbind(T,ER,MT )
ss <- cbind(ss,as.numeric(gsub(".*\\.","",row.names(ss))))
ss <- cbind(ss,comx)
ss <- ss[order(ss[,4]),]
ss <- cbind(ss,row.names(ss))
row.names(ss) <- ss[,4]


barplot(t(ss[match(id[,1], ss[,4]),1:3]),col=c("grey85","magenta","orange") ,border=NA ,ylab="Absolute Number",xlab="Single Cells (1:384)",main=paste("Plot of absolute reads",j,sep="\n"))
plot.new()
legend("left",c("Unique mRNA","ERCC","Mitochondrial RNA"),col=c("grey85","magenta","orange") , cex=1.5,pt.cex=3,bty="n",pch=15 )

pp <- cbind(U/T*100,ER/T*100,MT/T*100 )
pp <- cbind(pp,as.numeric(gsub(".*\\.","",row.names(pp))))
pp <- pp[order(pp[,4]),]
pp <- cbind(pp,row.names(pp))
row.names(pp) <- pp[,4]

barplot(t(pp[match(id[,1], pp[,4]),1:3]),col=c("dodgerblue","magenta","orange") ,border=NA,xlab="Plate barcodes", ylab="Percentage",main=paste("Proportional plot of the reads",j,sep="\n") )
plot.new()
legend("left",c("Unique mRNA","ERCC","Mitochondrial RNA"),col=c("dodgerblue","magenta","orange") , cex=1.5,pt.cex=3,bty="n",pch=15 )


## %UMI
q <- cbind(pp[match(id[,1], pp[,4]),],rep(16:1,each=24),rep(1:24,16))
col <- colorRampPalette(c("snow1","blue"))

q <- cbind(q,col(10)[cut(as.numeric(as.character(q[,1])),breaks = 10)])

plot(x=q[,7],y=q[,6],pch=19,cex=2,yaxt="n",col=ifelse(q[,1]==0|q[,1]=="NaN","red",q[,8]),xlab="",ylab="",main=paste("Percentage of mRNA molecules",j,sep="\n"),xaxt="n")
axis(1,at=c(1:24),labels=seq(1:24),cex.axis=0.5)
axis(2,at=c(16:1),labels=rep(LETTERS,length.out=16),cex.axis=0.5)

barplot(as.matrix(rep(1,11)),col=c("red",col(10)),border=NA,yaxt="n",width=2,xlim=c(0,20))
axis(2,at=c(.5,1.5:10.5),labels=c("0",gsub("\\,","\\-",levels(cut(as.numeric(as.character(q[,1])),breaks = 10)))))

## %ERCC
q <- cbind(pp[match(id[,1], pp[,4]),],rep(16:1,each=24),rep(1:24,16))
col <- colorRampPalette(c("snow1","blue"))

q <- cbind(q,col(10)[cut(as.numeric(as.character(q[,2])),breaks = 10)])

plot(x=q[,7],y=q[,6],pch=19,cex=2,yaxt="n",col=ifelse(q[,2]==0|q[,2]=="NaN","red",q[,8]),xlab="",ylab="",main=paste("Percentage of ERCC molecules",j,sep="\n"),xaxt="n")
axis(1,at=c(1:24),labels=seq(1:24),cex.axis=0.5)
axis(2,at=c(16:1),labels=rep(LETTERS,length.out=16),cex.axis=0.5)

barplot(as.matrix(rep(1,11)),col=c("red",col(10)),border=NA,yaxt="n",width=2,xlim=c(0,20))
axis(2,at=c(.5,1.5:10.5),labels=c("0",gsub("\\,","\\-",levels(cut(as.numeric(as.character(q[,2])),breaks = 10)))))

## %Mito
q <- cbind(pp[match(id[,1], pp[,4]),],rep(16:1,each=24),rep(1:24,16))
col <- colorRampPalette(c("snow1","blue"))

q <- cbind(q,col(10)[cut(as.numeric(as.character(q[,3])),breaks = 10)])

plot(x=q[,7],y=q[,6],pch=19,cex=2,yaxt="n",col=ifelse(q[,3]==0|q[,3]=="NaN","red",q[,8]),xlab="",ylab="",main=paste("Percentage of Mitochondrial mRNA molecules",j,sep="\n"),xaxt="n")
axis(1,at=c(1:24),labels=seq(1:24),cex.axis=0.5)
axis(2,at=c(16:1),labels=rep(LETTERS,length.out=16),cex.axis=0.5)

barplot(as.matrix(rep(1,11)),col=c("red",col(10)),border=NA,yaxt="n",width=2,xlim=c(0,20))
axis(2,at=c(.5,1.5:10.5),labels=c("0",gsub("\\,","\\-",levels(cut(as.numeric(as.character(q[,3])),breaks = 10)))))

##################################

## n-UMI
q <- cbind(ss[match(id[,1], ss[,4]),],rep(16:1,each=24),rep(1:24,16))
col <- colorRampPalette(c("snow1","blue"))

q <- cbind(q,col(10)[cut(as.numeric(as.character(q[,1])),breaks = 10)])

plot(x=q[,8],y=q[,7],pch=19,cex=2,yaxt="n",col=ifelse(q[,1]==0,"red",q[,9]),xlab="",ylab="",main=paste("Number of mRNA molecules",j,sep="\n"),xaxt="n")
axis(1,at=c(1:24),labels=seq(1:24),cex.axis=0.5)
axis(2,at=c(16:1),labels=rep(LETTERS,length.out=16),cex.axis=0.5)

barplot(as.matrix(rep(1,11)),col=c("red",col(10)),border=NA,yaxt="n",width=2,xlim=c(0,20))
axis(2,at=c(.5,1.5:10.5),labels=c("0",gsub("\\,","\\-",levels(cut(as.numeric(as.character(q[,1])),breaks = 10)))))

## n-ERCC
q <- cbind(ss[match(id[,1], ss[,4]),],rep(16:1,each=24),rep(1:24,16))
col <- colorRampPalette(c("snow1","blue"))

q <- cbind(q,col(10)[cut(as.numeric(as.character(q[,2])),breaks = 10)])

plot(x=q[,8],y=q[,7],pch=19,cex=2,yaxt="n",col=ifelse(q[,2]==0,"red",q[,9]),xlab="",ylab="",main=paste("Number of ERCC molecules",j,sep="\n"),xaxt="n")
axis(1,at=c(1:24),labels=seq(1:24),cex.axis=0.5)
axis(2,at=c(16:1),labels=rep(LETTERS,length.out=16),cex.axis=0.5)

barplot(as.matrix(rep(1,11)),col=c("red",col(10)),border=NA,yaxt="n",width=2,xlim=c(0,20))
axis(2,at=c(.5,1.5:10.5),labels=c("0",gsub("\\,","\\-",levels(cut(as.numeric(as.character(q[,2])),breaks = 10)))))

## n-Mito
q <- cbind(ss[match(id[,1], ss[,4]),],rep(16:1,each=24),rep(1:24,16))
col <- colorRampPalette(c("snow1","blue"))

q <- cbind(q,col(10)[cut(as.numeric(as.character(q[,3])),breaks = 10)])

plot(x=q[,8],y=q[,7],pch=19,cex=2,yaxt="n",col=ifelse(q[,3]==0,"red",q[,9]),xlab="",ylab="",main=paste("Number of Mitochondrial mRNA molecules",j,sep="\n"),xaxt="n")
axis(1,at=c(1:24),labels=seq(1:24),cex.axis=0.5)
axis(2,at=c(16:1),labels=rep(LETTERS,length.out=16),cex.axis=0.5)

barplot(as.matrix(rep(1,11)),col=c("red",col(10)),border=NA,yaxt="n",width=2,xlim=c(0,20))
axis(2,at=c(.5,1.5:10.5),labels=c("0",gsub("\\,","\\-",levels(cut(as.numeric(as.character(q[,3])),breaks = 10)))))

## gene coverage

q <- cbind(ss[match(id[,1], ss[,4]),],rep(16:1,each=24),rep(1:24,16))
col <- colorRampPalette(c("snow1","blue"))

q <- cbind(q,col(10)[cut(as.numeric(as.character(q[,5])),breaks = 10)])

plot(x=q[,8],y=q[,7],pch=19,cex=2,yaxt="n",col=ifelse(q[,5]==0,"red",q[,9]),xlab="",ylab="",main=paste("Gene coverage (n >= 1)",j,sep="\n"),xaxt="n")
axis(1,at=c(1:24),labels=seq(1:24),cex.axis=0.5)
axis(2,at=c(16:1),labels=rep(LETTERS,length.out=16),cex.axis=0.5)

barplot(as.matrix(rep(1,11)),col=c("red",col(10)),border=NA,yaxt="n",width=2,xlim=c(0,20))
axis(2,at=c(.5,1.5:10.5),labels=c("0",gsub("\\,","\\-",levels(cut(as.numeric(as.character(q[,5])),breaks = 10)))))

hist(comx,breaks=15,col="grey",xlab="Gene Coverage",main=paste("Gene Complexity Histogram",j,sep="\n"))

hist(ERcomx,breaks=15,col="grey40",xlab="ERCC Coverage",main=paste("ERCC Complexity Histogram",j,sep="\n"))

while(!par('page')) plot.new()

}

dev.off()

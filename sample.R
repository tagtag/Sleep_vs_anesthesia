#--- short read mapping ---------------
#all short read mapping was done the following command (as an exmple, DRR231827 is considered)
# bowtie2-2.3.5.1-linux-x86_64/bowtie2 -x ~/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome -1 DRR231827_1.fastq.bz2 -2 DRR231827_2.fastq.bz2 -p 12 -S DRR231827.sam
# samtools view  -Sb DRR231827.sam | samtools sort - -@ 12 -o DRR231827_sorted.bam 
# samtools index DRR231827_sorted.bam
# htseq-count -f bam DRR231827_sorted.bam ~/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf > DRR231827.txt 
# all corresoponding .txt files are in ./data folder
# 

#  load GSE69079_series_matrix.txt.gz and GPL1261-56135.txt.gz dowloaded from GEO
x <- read.csv("GSE69079_series_matrix.txt.gz",sep="\t",comment.char="!")
y <- read.csv("GPL1261-56135.txt.gz",sep="\t",comment.char="#")
# load Bioconductor package org.Mm.eg.db
require("org.Mm.eg.db")
# prepare a file
x1 <- read.csv("DRR231827.txt",sep="\t",header=F)
gene_info <- select(org.Mm.eg.db,
                    keys = x1[,1],
                    columns = c("SYMBOL", "ENTREZID"),
                    keytype = "SYMBOL")

index <- match(y$ID[match(gene_info[,2],y$ENTREZ_GENE_ID)],x[,1])
# load all .txt files and combine them into one data frame
files <- list.files(pattern = "^DRR.*\\.txt$")
x1_all <- NULL
for (i in seq_along(files))
{
  cat(i," ")
  x1 <- read.csv(files[i],sep="\t",header=F)
  if (i==1) {
    x1_all <- x1
  } else {
      x1_all <- cbind(x1_all,x1[,2])
    }
}
# construct tensor 
Z1 <- array(NA,c(dim(x1_all)[1],6,2,2))
Z1[,,1,1] <- data.matrix(x1_all[,2:7])
Z1[,,1,2] <-  data.matrix(x1_all[,14:19])
Z1[,,2,1] <-  data.matrix(x1_all[,rep(20:21,each=3)])
Z1[,,2,2] <-  data.matrix(x1_all[,rep(24:25,each=3)])

Z <- array(NA,c(dim(x)[1],3,4,2))
Z[,,1,1] <- data.matrix(x[,1+c(1,3,5)])
Z[,,1,2] <- data.matrix(x[,1+c(2,4,6)])
Z[,,2,1] <- data.matrix(x[,7+c(1,3,5)])
Z[,,2,2] <- data.matrix(x[,7+c(2,4,6)])
Z[,,3,1] <- data.matrix(x[,13+c(1,3,5)])
Z[,,3,2] <- data.matrix(x[,13+c(2,4,6)])
Z[,,4,1] <- data.matrix(x[,19+c(1,3,5)])
Z[,,4,2] <- data.matrix(x[,19+c(2,4,6)])

Z1 <- Z1[!is.na(index),,,]
genes <- x1[!is.na(index),1]
index <- index[!is.na(index)]

ZZ <- array(NA,c(dim(Z1)[1:3],dim(Z)[2:4]))
for(i in seq_len(dim(ZZ)[1]))
{
  cat(i," ")
  ZZ_tmp <- outer(Z1[i,,,],Z[index[i],,,])
 ZZ[i,,,,,1] <-ZZ_tmp[,,1,,,1]
 ZZ[i,,,,,2] <-ZZ_tmp[,,2,,,2]
}
ZZ <- apply(ZZ,2:6,scale)
# perform tensodr decomposition
require(rTensor)
HOSVD <- hosvd(as.tensor(ZZ),c(10,6,2,3,4,2))
# select genes using k0<-3 (l_1=3) or k0<-4 (l_1=4)
k0<-3
#k0<-4
th <- function(sd){
  P2<- pchisq(((HOSVD$U[[1]][,k0]-mean(HOSVD$U[[1]][,k0]))/sd)^2,1,lower.tail=F)
  hc<- hist(1-P2,breaks=100,plot=F)
  return(sd(hc$count[1:sum(hc$mid<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.001,th)$par 

aa <- seq(0.5*sd,2*sd,by=0.05*sd)
bb<-apply(matrix(seq(0.5*sd,2*sd,by=0.05*sd),ncol=1),1,th)
pdf(file="hist.pdf",width=10,height=5) #k0<-3
#pdf(file="hist_2.pdf",width=10,height=5) #k0<-4
par(mfrow=c(1,2))
plot(aa,bb,xlab="sigma_l",ylab="sigma_h",type="o",cex.lab=2,cex.axis=2)
arrows(sd,max(bb),sd,min(bb),col=2)
hist(1-P1,breaks=100,xlab="1-Pi",cex.lab=2,cex.axis=2)
par(mfrow=c(1,1))
dev.off()

genes_plus <- genes[p.adjust(P1,"BH")<0.01 & HOSVD$U[[1]][,k0]>0 ] #k0<-3
#genes_minus <- genes[p.adjust(P1,"BH")<0.01 & HOSVD$U[[1]][,k0]<0 ] #k0<-4
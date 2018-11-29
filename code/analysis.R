###################################################################################################
# SETTING UP THE ENVIRONEMNT
###################################################################################################

library(readr)
library(corrplot)
library(boot)
library(ggplot2)
library(reshape)
library(here)
library(ppcor)

###################################################################################################
# HELPER FUNCTIONS
###################################################################################################

corrFun <- function(data,indices){
  d <- data[indices,] # alows boot to select samples
  return(cor.test(d[,1],d[,2], method = "spearman")$estimate[[1]])
}

###################################################################################################
# GENES IN EACH PANEL
###################################################################################################

decipher_genes <- c("LASP1","IQGAP3","NFIB","S1PR4","THBS2","ANO7","PCDH7","MYBPC1","EPPK1","TSBP",
                    "PBX1","NUSAP1","ZWILCH","UBE2C","CAMK2N1","RABGAP1","PCAT-32","PCAT-80","TNFRSF19","C6orf10")

prolaris_genes <- c("FOXM1","CDC20","CDKN3","CDC2","KIF11","KIAA0101","NUSAP1","CENPF","ASPM",
                    "BUB1B","RRM2","DLGAP5","BIRC5","KIF20A","PLK1","TOP2A","TK1","PBK","ASF1B",
                    "C18orf24","RAD54L","PTTG1","CDCA3","MCM10","PRC1","DTL","CEP55","RAD51",
                    "CENPM","CDCA8","ORC6L","SKA1","ORC6","CDK1") 

oncotypedx_genes <- c("AZGP1","KLK2","SRD5A2","FAM13C","FLNC","GSN","TPM2","GSTM2","TPX2","BGN",
                      "COL1A1","SFRP4")

###################################################################################################
# DATA IMPORT
###################################################################################################


dat <- read_delim(here("data/nano.txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
dat <- as.data.frame(dat)

dat$CAPRASn <- ifelse(dat$CAPRASgroup=='LOW',0,
                      ifelse(dat$CAPRASgroup=="HIGH",2,1))

###################################################################################################
# DIFFERNETIAL XPRESSION - EAM VS AAM
###################################################################################################

# calculate p value for mann whitney between races for each gene
raceMannAll <- apply(dat[,colnames(dat) %in% c(prolaris_genes,oncotypedx_genes,decipher_genes)],
                     2,function(x) wilcox.test(x~dat$Race)$p.value) 
raceMannAll <- as.data.frame(cbind(gene=colnames(dat)[2:61],
                                   raceMannAll))
raceMannAll$raceMannAll <- as.numeric(as.character(raceMannAll$raceMannAll))

# annoate labels
raceMannAll$label <- ifelse(raceMannAll$raceMannAll<=0.001,"***",
                            ifelse(raceMannAll$raceMannAll<=0.01,"**",
                                   ifelse(raceMannAll$raceMannAll<=0.05,"*",NA)))

# calculate fold change
raceMannAll$fold <- NA

for(i in 1:length(raceMannAll$gene)){
  raceMannAll$fold[i] <- round((median(dat[dat$Race=="EAM",as.character(raceMannAll$gene[i])])/median(dat[dat$Race=="AAM",as.character(raceMannAll$gene[i])])),2)
}

raceMannAll$fold <- round(log(raceMannAll$fold,2),2)
raceMannAll <- raceMannAll[!is.na(raceMannAll$label),]

melted_data <- dat[,1:62]
melted_data <- melt(melted_data, id.vars = c("ID","Race"))

### Figure 1A
ggplot(melted_data[melted_data$variable %in% prolaris_genes,], aes(x=variable, y=value, fill=as.factor(Race)))+
  scale_fill_manual(values=rep("white",length(prolaris_genes)))+
  geom_jitter(position=position_jitterdodge(jitter.width=.15, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.25)+
  geom_boxplot(outlier.shape=NA,alpha=0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="Gene",y="Expression",title="Prolaris") +
  annotate("text",label=raceMannAll[raceMannAll$gene %in% prolaris_genes,]$label,x=raceMannAll[raceMannAll$gene %in% prolaris_genes,]$gene,y=9.5) + 
  annotate("text",label=raceMannAll[raceMannAll$gene %in% prolaris_genes,]$fold,x=raceMannAll[raceMannAll$gene %in% prolaris_genes,]$gene,y=0, cex=3.5) + 
  scale_color_manual(values=c("deepskyblue", "tomato"))

### Figure 1B
ggplot(melted_data[melted_data$variable %in% oncotypedx_genes,], aes(x=variable, y=value, fill=as.factor(Race)))+
  scale_fill_manual(values=rep("white",length(oncotypedx_genes)))+
  geom_jitter(position=position_jitterdodge(jitter.width=.15, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.25)+
  geom_boxplot(outlier.shape=NA, alpha=0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15),
        legend.text=element_text(size=15), legend.title=element_text(size=15)) +
  guides(fill=FALSE, colour=guide_legend(title="Race"), colour=FALSE) +
  labs(x="Gene",y="Expression",title="Oncotype DX") +
  annotate("text",label=raceMannAll[raceMannAll$gene %in% oncotypedx_genes,]$label,x=raceMannAll[raceMannAll$gene %in% oncotypedx_genes,]$gene,y=17) +
  annotate("text",label=raceMannAll[raceMannAll$gene %in% oncotypedx_genes,]$fold,x=raceMannAll[raceMannAll$gene %in% oncotypedx_genes,]$gene,y=0, cex=3.5) + 
  scale_color_manual(values=c("deepskyblue", "tomato"), labels=c("AAM","EAM"))

### Figure 1C
ggplot(melted_data[melted_data$variable %in% decipher_genes,], aes(x=variable, y=value, fill=as.factor(Race)))+
  scale_fill_manual(values=rep("white",length(decipher_genes)))+
  geom_jitter(position=position_jitterdodge(jitter.width=.15, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.25)+
  geom_boxplot(outlier.shape=NA, alpha=0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, face = "italic", color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="Gene",y="Expression",title="Decipher") +
  annotate("text",label=raceMannAll[raceMannAll$gene %in% decipher_genes,]$label,x=raceMannAll[raceMannAll$gene %in% decipher_genes,]$gene,y=15) + 
  annotate("text",label=raceMannAll[raceMannAll$gene %in% decipher_genes,]$fold,x=raceMannAll[raceMannAll$gene %in% decipher_genes,]$gene,y=0, 
           cex=3.5) + 
  scale_color_manual(values=c("deepskyblue", "tomato"))

###################################################################################################
# INTER-GENE CORRELATION (SEPERATELY IN EAM AND AAM) AND ICC
###################################################################################################

# Prolaris

## Calculate spearmans correlations rounded to one decimal place
corrAAM <- round(cor(dat[dat$Race=="AAM",colnames(dat) %in% prolaris_genes],method="spearman"),1)
corrEAM <- round(cor(dat[dat$Race=="EAM",colnames(dat) %in% prolaris_genes],method="spearman"),1)

### Figure 2A1
corrplot(corrEAM[c("KIAA0101","RAD54L","SKA1","CDCA3","CENPM","RRM2","KIF11","PTTG1","KIF20A",
                   "FOXM1","TK1","ASF1B","MCM10","PBK","ASPM","CDC20","BUB1B","NUSAP1","CEP55",
                   "BIRC5","CENPF","PLK1","PRC1","RAD51","ORC6","DTL","DLGAP5","CDK1",
                   "CDKN3","TOP2A"),
                 c("KIAA0101","RAD54L","SKA1","CDCA3","CENPM","RRM2","KIF11","PTTG1","KIF20A",
                   "FOXM1","TK1","ASF1B","MCM10","PBK","ASPM","CDC20","BUB1B","NUSAP1","CEP55",
                   "BIRC5","CENPF","PLK1","PRC1","RAD51","ORC6","DTL","DLGAP5","CDK1",
                   "CDKN3","TOP2A")],
         cl.lim = c(-0.5,1), method="color", col= colorRampPalette(c("blue","white", "red"))(20),
         tl.cex=1.3, tl.col = "black",cl.cex=0.9, cl.ratio = 0.1,addgrid.col = "lightgrey",
         family="sans",cl.pos="n", font=3)

### Figure 2A2
corrplot(corrAAM[c("KIAA0101","RAD54L","SKA1","CDCA3","CENPM","RRM2","KIF11","PTTG1","KIF20A",
                   "FOXM1","TK1","ASF1B","MCM10","PBK","ASPM","CDC20","BUB1B","NUSAP1","CEP55",
                   "BIRC5","CENPF","PLK1","PRC1","RAD51","ORC6","DTL","DLGAP5","CDK1",
                   "CDKN3","TOP2A"),
                 c("KIAA0101","RAD54L","SKA1","CDCA3","CENPM","RRM2","KIF11","PTTG1","KIF20A",
                   "FOXM1","TK1","ASF1B","MCM10","PBK","ASPM","CDC20","BUB1B","NUSAP1","CEP55",
                   "BIRC5","CENPF","PLK1","PRC1","RAD51","ORC6","DTL","DLGAP5","CDK1",
                   "CDKN3","TOP2A")],
         cl.lim = c(-0.5,1), method="color", col= colorRampPalette(c("blue","white", "red"))(20),
         tl.cex=1.3, tl.col = "black",cl.cex=0.9, cl.ratio = 0.1,addgrid.col = "lightgrey",
         family="sans",cl.pos="n", font=3)

## Calculate ICC
## grab upper triangle on correlation values - recalculate corr without rounding 
corrAAM <- cor(dat[dat$Race=="AAM",colnames(dat) %in% prolaris_genes],method="spearman")
corrEAM <- cor(dat[dat$Race=="EAM",colnames(dat) %in% prolaris_genes],method="spearman")

vectorAlpha <- corrAAM[upper.tri(corrAAM)]
vectorBeta <- corrEAM[upper.tri(corrEAM)]

## calculate spearmans correlation of correlations
cor.test(vectorAlpha,vectorBeta, method = "spearman")

## make data set for bootstrapping and plotting 
qdat <- as.data.frame(cbind(vectorAlpha,vectorBeta))

## bootstrap for CI
boot_results <- boot(data=qdat, statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

## Calculate ICC - CAPRAS weighted
corrAAM2 <- c()
for ( i in 1:ncol(dat[,colnames(dat) %in% prolaris_genes])){
  for ( j in 1:ncol(dat[,colnames(dat) %in% prolaris_genes])){
    if ( i == j ) { corrAAM2 <- c(corrAAM2,1)}
    else{
      corrAAM2 <- c(corrAAM2,pcor.test(dat[dat$Race=="AAM",colnames(dat) %in% prolaris_genes][,i],
                                       dat[dat$Race=="AAM",colnames(dat) %in% prolaris_genes][,j],
                                       dat[dat$Race=="AAM",]$CAPRASn, method="spearman")$estimate)
    }
  }
}

corrAAM2 <- matrix(corrAAM2, ncol=ncol(dat[,colnames(dat) %in% prolaris_genes]))

corrEAM2 <- c()
for ( i in 1:ncol(dat[,colnames(dat) %in% prolaris_genes])){
  for ( j in 1:ncol(dat[,colnames(dat) %in% prolaris_genes])){
    if ( i == j ) { corrEAM2 <- c(corrEAM2,1)}
    else{
      corrEAM2 <- c(corrEAM2,pcor.test(dat[dat$Race=="EAM",colnames(dat) %in% prolaris_genes][,i],
                                       dat[dat$Race=="EAM",colnames(dat) %in% prolaris_genes][,j],
                                       dat[dat$Race=="EAM",]$CAPRASn, method="spearman")$estimate)
    }
  }
}

corrEAM2 <- matrix(corrEAM2, ncol=ncol(dat[,colnames(dat) %in% prolaris_genes]))

# Step 2 : make a vector
vectorAlpha <- corrAAM2[upper.tri(corrAAM2)]
vectorBeta <- corrEAM2[upper.tri(corrEAM2)]

# step 3: corr of corr
cor.test(vectorAlpha,vectorBeta, method = "spearman")

qdat <- as.data.frame(cbind(vectorAlpha,vectorBeta))

# bootstrap for CI
boot_results <- boot(data=qdat, statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")


### Figure 2A3
ggplot(qdat,aes(x=vectorAlpha, y=vectorBeta))+
  geom_point(size=3) +
  annotate("text",label="CORR = 0.62 (0.56 , 0.69)",x=-0.5,y=0.8,cex=5) +
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  xlim(-1,1) + ylim(-1,1) +
  geom_smooth(se=FALSE, colour="purple", method="glm",fullrange=TRUE)+
  labs(x="AAM",y="EAM", title="Prolaris")

###################################################################################################
# Oncotype DX

## Calculate spearmans correlations rounded to one decimal place
corrAAM <- round(cor(dat[dat$Race=="AAM",colnames(dat) %in% oncotypedx_genes],method="spearman"),1)
corrEAM <- round(cor(dat[dat$Race=="EAM",colnames(dat) %in% oncotypedx_genes],method="spearman"),1)

### Figure 2B1
corrplot(corrEAM[c("TPX2","SFRP4","BGN","COL1A1","AZGP1","FAM13C","KLK2","SRD5A2","GSN","GSTM2","FLNC","TPM2"),
                 c("TPX2","SFRP4","BGN","COL1A1","AZGP1","FAM13C","KLK2","SRD5A2","GSN","GSTM2","FLNC","TPM2")],
         cl.lim = c(-0.5,1), method="color", col= colorRampPalette(c("blue","white", "red"))(20),
         tl.cex=1.5, tl.col = "black",cl.cex=0.6, cl.ratio = 0.1,addgrid.col = "lightgrey",
         family="sans", font=3)

### Figure 2B2
corrplot(corrAAM[c("TPX2","SFRP4","BGN","COL1A1","AZGP1","FAM13C","KLK2","SRD5A2","GSN","GSTM2","FLNC","TPM2"),
                 c("TPX2","SFRP4","BGN","COL1A1","AZGP1","FAM13C","KLK2","SRD5A2","GSN","GSTM2","FLNC","TPM2")],
         cl.lim = c(-0.5,1), method="color", col= colorRampPalette(c("blue","white", "red"))(20),
         tl.cex=1.5, tl.col = "black",cl.cex=0.6, cl.ratio = 0.1,addgrid.col = "lightgrey",
         family="sans", font=3)

## Calculate ICC
## grab upper triangle on correlation values - recalculate corr without rounding 
corrAAM <- cor(dat[dat$Race=="AAM",colnames(dat) %in% oncotypedx_genes],method="spearman")
corrEAM <- cor(dat[dat$Race=="EAM",colnames(dat) %in% oncotypedx_genes],method="spearman")

vectorAlpha <- corrAAM[upper.tri(corrAAM)]
vectorBeta <- corrEAM[upper.tri(corrEAM)]

## calculate spearmans correlation of correlations
cor.test(vectorAlpha,vectorBeta, method = "spearman")

## make data set for bootstrapping and plotting 
qdat <- as.data.frame(cbind(vectorAlpha,vectorBeta))

## bootstrap for CI
boot_results <- boot(data=qdat, statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

## Calculate ICC - CAPRAS weighted
corrAAM2 <- c()
for ( i in 1:ncol(dat[,colnames(dat) %in% oncotypedx_genes])){
  for ( j in 1:ncol(dat[,colnames(dat) %in% oncotypedx_genes])){
    if ( i == j ) { corrAAM2 <- c(corrAAM2,1)}
    else{
      corrAAM2 <- c(corrAAM2,pcor.test(dat[dat$Race=="AAM",colnames(dat) %in% oncotypedx_genes][,i],
                                       dat[dat$Race=="AAM",colnames(dat) %in% oncotypedx_genes][,j],
                                       dat[dat$Race=="AAM",]$CAPRASn, method="spearman")$estimate)
    }
  }
}

corrAAM2 <- matrix(corrAAM2, ncol=ncol(dat[,colnames(dat) %in% oncotypedx_genes]))

corrEAM2 <- c()
for ( i in 1:ncol(dat[,colnames(dat) %in% oncotypedx_genes])){
  for ( j in 1:ncol(dat[,colnames(dat) %in% oncotypedx_genes])){
    if ( i == j ) { corrEAM2 <- c(corrEAM2,1)}
    else{
      corrEAM2 <- c(corrEAM2,pcor.test(dat[dat$Race=="EAM",colnames(dat) %in% oncotypedx_genes][,i],
                                       dat[dat$Race=="EAM",colnames(dat) %in% oncotypedx_genes][,j],
                                       dat[dat$Race=="EAM",]$CAPRASn, method="spearman")$estimate)
    }
  }
}

corrEAM2 <- matrix(corrEAM2, ncol=ncol(dat[,colnames(dat) %in% oncotypedx_genes]))

# Step 2 : make a vector
vectorAlpha <- corrAAM2[upper.tri(corrAAM2)]
vectorBeta <- corrEAM2[upper.tri(corrEAM2)]

# step 3: corr of corr
cor.test(vectorAlpha,vectorBeta, method = "spearman")

qdat <- as.data.frame(cbind(vectorAlpha,vectorBeta))

# bootstrap for CI
boot_results <- boot(data=qdat, statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

### Figure 2B3
ggplot(qdat,aes(x=vectorAlpha, y=vectorBeta))+
  geom_point(size=3) +
  annotate("text",label="CORR = 0.87 (0.81 , 0.97)",x=-0.5,y=0.8,cex=5) +
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  xlim(-1,1) + ylim(-1,1) +
  geom_smooth(se=FALSE, colour="purple", method="glm",fullrange=TRUE)+
  labs(x="AAM",y="EAM", title="Oncotype")

###################################################################################################
# Decipher

## Calculate spearmans correlations rounded to one decimal place
corrAAM <- round(cor(dat[dat$Race=="AAM",colnames(dat) %in% decipher_genes],method="spearman"),1)
corrEAM <- round(cor(dat[dat$Race=="EAM",colnames(dat) %in% decipher_genes],method="spearman"),1)

### Figure 2C1
corrplot(corrEAM[c("PBX1","PCDH7","EPPK1","CAMK2N1","THBS2","ZWILCH","S1PR4","C6orf10","NUSAP1",
                   "UBE2C","IQGAP3","PCAT-32","PCAT-80","TNFRSF19","ANO7","LASP1","RABGAP1",
                   "MYBPC1","NFIB"),
                 c("PBX1","PCDH7","EPPK1","CAMK2N1","THBS2","ZWILCH","S1PR4","C6orf10","NUSAP1",
                   "UBE2C","IQGAP3","PCAT-32","PCAT-80","TNFRSF19","ANO7","LASP1","RABGAP1",
                   "MYBPC1","NFIB")],
         cl.lim = c(-0.5,1), method="color", col= colorRampPalette(c("blue","white", "red"))(20),
         tl.cex=1.3, tl.col = "black",cl.cex=0.6, cl.ratio = 0.1,addgrid.col = "lightgrey",
         family="sans",cl.pos="n",font=3)

### Figure 2C2
corrplot(corrAAM[c("PBX1","PCDH7","EPPK1","CAMK2N1","THBS2","ZWILCH","S1PR4","C6orf10","NUSAP1",
                   "UBE2C","IQGAP3","PCAT-32","PCAT-80","TNFRSF19","ANO7","LASP1","RABGAP1",
                   "MYBPC1","NFIB"),
                 c("PBX1","PCDH7","EPPK1","CAMK2N1","THBS2","ZWILCH","S1PR4","C6orf10","NUSAP1",
                   "UBE2C","IQGAP3","PCAT-32","PCAT-80","TNFRSF19","ANO7","LASP1","RABGAP1",
                   "MYBPC1","NFIB")],
         cl.lim = c(-0.5,1), method="color", col= colorRampPalette(c("blue","white", "red"))(20),
         tl.cex=1.3, tl.col = "black",cl.cex=0.6, cl.ratio = 0.1,addgrid.col = "lightgrey",
         family="sans",cl.pos="n",font=3)

## Calculate ICC
## grab upper triangle on correlation values - recalculate corr without rounding 
corrAAM <- cor(dat[dat$Race=="AAM",colnames(dat) %in% decipher_genes],method="spearman")
corrEAM <- cor(dat[dat$Race=="EAM",colnames(dat) %in% decipher_genes],method="spearman")

vectorAlpha <- corrAAM[upper.tri(corrAAM)]
vectorBeta <- corrEAM[upper.tri(corrEAM)]

## calculate spearmans correlation of correlations
cor.test(vectorAlpha,vectorBeta, method = "spearman")

## make data set for bootstrapping and plotting 
qdat <- as.data.frame(cbind(vectorAlpha,vectorBeta))

## bootstrap for CI
boot_results <- boot(data=qdat, statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

## Calculate ICC - CAPRAS weighted
corrAAM2 <- c()
for ( i in 1:ncol(dat[,colnames(dat) %in% decipher_genes])){
  for ( j in 1:ncol(dat[,colnames(dat) %in% decipher_genes])){
    if ( i == j ) { corrAAM2 <- c(corrAAM2,1)}
    else{
      corrAAM2 <- c(corrAAM2,pcor.test(dat[dat$Race=="AAM",colnames(dat) %in% decipher_genes][,i],
                                       dat[dat$Race=="AAM",colnames(dat) %in% decipher_genes][,j],
                                       dat[dat$Race=="AAM",]$CAPRASn, method="spearman")$estimate)
    }
  }
}

corrAAM2 <- matrix(corrAAM2, ncol=ncol(dat[,colnames(dat) %in% decipher_genes]))

corrEAM2 <- c()
for ( i in 1:ncol(dat[,colnames(dat) %in% decipher_genes])){
  for ( j in 1:ncol(dat[,colnames(dat) %in% decipher_genes])){
    if ( i == j ) { corrEAM2 <- c(corrEAM2,1)}
    else{
      corrEAM2 <- c(corrEAM2,pcor.test(dat[dat$Race=="EAM",colnames(dat) %in% decipher_genes][,i],
                                       dat[dat$Race=="EAM",colnames(dat) %in% decipher_genes][,j],
                                       dat[dat$Race=="EAM",]$CAPRASn, method="spearman")$estimate)
    }
  }
}

corrEAM2 <- matrix(corrEAM2, ncol=ncol(dat[,colnames(dat) %in% decipher_genes]))

# Step 2 : make a vector
vectorAlpha <- corrAAM2[upper.tri(corrAAM2)]
vectorBeta <- corrEAM2[upper.tri(corrEAM2)]

# step 3: corr of corr
cor.test(vectorAlpha,vectorBeta, method = "spearman")

qdat <- as.data.frame(cbind(vectorAlpha,vectorBeta))

# bootstrap for CI
boot_results <- boot(data=qdat, statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

### Figure 2C3
ggplot(qdat,aes(x=vectorAlpha, y=vectorBeta))+
  geom_point(size=3) +
  annotate("text",label="CORR = 0.73 (0.66 , 0.82)",x=-0.5,y=0.8,cex=5) +
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  xlim(-1,1) + ylim(-1,1) +
  geom_smooth(se=FALSE, colour="purple", method="glm",fullrange=TRUE)+
  labs(x="AAM",y="EAM", title="Decipher")

###################################################################################################
# RISK SCORES - EAM VS AAM
###################################################################################################
# Prolaris

## calculate risk score
prolarisScore <- apply(dat[colnames(dat) %in% prolaris_genes],1,function(x) sum(x)) 

risk_prolaris <- as.data.frame(cbind(id=dat$ID,
                                     score=prolarisScore,
                                     Race=dat$Race,
                                     capra=dat$CAPRASgroup,
                                     nccn=dat$nccn))

risk_prolaris$score <- as.numeric(as.character(risk_prolaris$score))

wilcox.test(risk_prolaris$score~risk_prolaris$Race)

### Figure 3B 
ggplot(risk_prolaris, aes(x=capra, y=score, fill=Race))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=rep("white",3))+
  geom_jitter(position=position_jitterdodge(jitter.width=0.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.7,size=3)+
  scale_x_discrete(limits = c("LOW","INTERMEDIATE","HIGH")) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="CAPRA-S",y="Risk Score",title="Prolaris") +
  scale_color_manual(values=c("deepskyblue", "tomato"))

### Figure 3B2 
ggplot(risk_prolaris, aes(x=nccn, y=score, fill=Race))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=rep("white",3))+
  geom_jitter(position=position_jitterdodge(jitter.width=0.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.7,size=3)+
  scale_x_discrete(limits = c("Low","Intermediate","High")) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="NCCN",y="Risk Score",title="Prolaris") +
  scale_color_manual(values=c("deepskyblue", "tomato"))

###################################################################################################
# oncotype dx

oncotypeScore <- dat$BGN + dat$COL1A1 + dat$SFRP4 - dat$FLNC + dat$GSN - dat$GSTM2 - dat$TPM2 -
  dat$AZGP1 - dat$FAM13C - dat$KLK2 - dat$SRD5A2 + dat$TPX2

risk_oncotype <- as.data.frame(cbind(id=dat$ID,
                                     score=oncotypeScore,
                                     Race=dat$Race,
                                     capra=dat$CAPRASgroup,
                                     nccn=dat$nccn))
risk_oncotype$score <- as.numeric(as.character(risk_oncotype$score))

wilcox.test(risk_oncotype$score~risk_oncotype$Race)

## Figure 3C
ggplot(risk_oncotype, aes(x=capra, y=score, fill=as.factor(Race)))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=rep("white",3))+
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.7,size=3)+
  scale_x_discrete(limits = c("LOW","INTERMEDIATE","HIGH")) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="CAPRA-S",y="Risk Score",title="Oncotype DX") +
  scale_color_manual(values=c("deepskyblue", "tomato"))

## Figure 3C2
ggplot(risk_oncotype, aes(x=nccn, y=score, fill=as.factor(Race)))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=rep("white",3))+
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.7,size=3)+
  scale_x_discrete(limits = c("Low","Intermediate","High")) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="NCCN",y="Risk Score",title="Oncotype DX") +
  scale_color_manual(values=c("deepskyblue", "tomato"))


###################################################################################################
# Decipher

decipherScore <- dat$LASP1 + dat$IQGAP3 + dat$NFIB -
  dat$S1PR4 + dat$THBS2 - dat$ANO7 - 
  dat$PCDH7 - dat$MYBPC1 + dat$EPPK1 -
  dat$C6orf10 + dat$PBX1 + dat$UBE2C + 
  dat$NUSAP1 + dat$ZWILCH + dat$CAMK2N1 +
  dat$RABGAP1 - dat$`PCAT-32` - dat$`PCAT-80` - dat$TNFRSF19


risk_decipher <- as.data.frame(cbind(id=dat$ID,
                                     score=decipherScore,
                                     Race=dat$Race,
                                     capra=dat$CAPRASgroup,
                                     nccn=dat$nccn))
risk_decipher$score <- as.numeric(as.character(risk_decipher$score))

wilcox.test(risk_decipher$score~risk_decipher$Race)

### Figure 3D
ggplot(risk_decipher, aes(x=capra, y=score, fill=as.factor(Race)))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=rep("white",3))+
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.7,size=3)+
  scale_x_discrete(limits = c("LOW","INTERMEDIATE","HIGH")) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="CAPRA-S",y="Risk Score",title="Decipher") +
  scale_color_manual(values=c("deepskyblue", "tomato"))

### Figure 3D2
ggplot(risk_decipher, aes(x=nccn, y=score, fill=as.factor(Race)))+
  geom_boxplot(outlier.shape=NA)+
  scale_fill_manual(values=rep("white",3))+
  geom_jitter(position=position_jitterdodge(jitter.width=.5, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.7,size=3)+
  scale_x_discrete(limits = c("Low","Intermediate","High")) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  guides(fill=FALSE, colour=FALSE) +
  labs(x="NCCN",y="Risk Score",title="Decipher") +
  scale_color_manual(values=c("deepskyblue", "tomato"))

###################################################################################################

### Figure 3A

tmp <- as.data.frame(cbind(Race=c("AAM","AAM","AAM","EAM","EAM","EAM"),
                           capra=c("LOW","INTERMEDIATE","HIGH","LOW","INTERMEDIATE","HIGH"),
                           perc=c(65,29,5,60,30,10)))

ggplot(tmp, aes(fill=Race, x=capra)) +
  geom_col(aes(y=as.numeric(as.character(perc))),position = "dodge") + 
  scale_x_discrete(limits = c("LOW","INTERMEDIATE","HIGH")) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        text = element_text(family="Calibri"),
        axis.title = element_text(face="bold", size=15)) +
  labs(x="CAPRA-S",y="%") +
  guides(fill=guide_legend(title="Race")) +
  scale_fill_manual(values=c("deepskyblue", "tomato"))

### Figure 3A2

tmp <- as.data.frame(cbind(Race=c("AAM","AAM","AAM","EAM","EAM","EAM"),
                           nccn=c("Low","Intermediate","High","Low","Intermediate","High"),
                           perc=c(4,77,19,9,54,37)))

ggplot(tmp, aes(fill=Race, x=nccn)) +
  geom_col(aes(y=as.numeric(as.character(perc))),position = "dodge") + 
  scale_x_discrete(limits = c("Low","Intermediate","High")) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        text = element_text(family="Calibri"),
        axis.title = element_text(face="bold", size=15)) +
  labs(x="NCCN",y="%") +
  guides(fill=guide_legend(title="Race")) +
  scale_fill_manual(values=c("deepskyblue", "tomato"))

###################################################################################################
# risk score correlations

# stratify risk scores by race
tmpE <- as.data.frame(cbind(decipher=risk_decipher[risk_decipher$Race=="EAM",]$score,
                            prolaris=risk_prolaris[risk_prolaris$Race=="EAM",]$score,
                            oncotype=risk_oncotype[risk_oncotype$Race=="EAM",]$score))

tmpA <- as.data.frame(cbind(decipher=risk_decipher[risk_decipher$Race=="AAM",]$score,
                            prolaris=risk_prolaris[risk_prolaris$Race=="AAM",]$score,
                            oncotype=risk_oncotype[risk_oncotype$Race=="AAM",]$score))

# decipher vs prolaris
cor.test(tmpE$decipher,tmpE$prolaris, method = "spearman")
cor.test(tmpA$decipher,tmpA$prolaris, method = "spearman")

# bootstrap for CI
boot_results <- boot(data=tmpA[,1:2], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")
boot_results <- boot(data=tmpE[,1:2], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

# Figure 4C
ggplot(tmpE,aes(x=decipher,y=prolaris))+
  geom_point(colour="tomato",size=3) +
  geom_point(data=tmpA,aes(x=decipher,y=prolaris),colour="deepskyblue",size=3) +
  annotate("text",label="EAM = -0.15 (-0.27 , -0.03)",x=45,y=110,cex=5) +
  annotate("text",label="AAM = 0.35 (0.17 , 0.55)",x=45,y=105,cex=5) +
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  geom_smooth(se=FALSE, colour="tomato", method="glm",fullrange=TRUE)+
  geom_smooth(data=tmpA,aes(x=decipher,y=prolaris),se=FALSE, colour="deepskyblue", method="glm",fullrange=TRUE)+
  labs(x="Decipher",y="Prolaris")

# decipher vs oncotype
cor.test(tmpE$decipher,tmpE$oncotype, method = "spearman")
cor.test(tmpA$decipher,tmpA$oncotype, method = "spearman")

# bootstrap for CI
boot_results <- boot(data=tmpA[,c(1,3)], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")
boot_results <- boot(data=tmpE[,c(1,3)], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

# Fifure 4B
ggplot(tmpE,aes(x=decipher,y=oncotype))+
  geom_point(colour="tomato",size=3) +
  geom_point(data=tmpA,aes(x=decipher,y=oncotype),colour="deepskyblue",size=3) +
  annotate("text",label="EAM = 0.39 (0.29 , 0.52)",x=27,y=-15,cex=5) +
  annotate("text",label="AAM = 0.49 (0.37 , 0.65)",x=27,y=-17,cex=5) +
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  geom_smooth(se=FALSE, colour="tomato", method="glm",fullrange=TRUE)+
  geom_smooth(data=tmpA,aes(x=decipher,y=oncotype),se=FALSE, colour="deepskyblue", method="glm",fullrange=TRUE)+
  labs(x="Decipher",y="Oncotype")



# prolaris vs oncotype
cor.test(tmpE$prolaris,tmpE$oncotype, method = "spearman")
cor.test(tmpA$prolaris,tmpA$oncotype, method = "spearman")

boot_results <- boot(data=tmpA[,2:3], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")
boot_results <- boot(data=tmpE[,2:3], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

ggplot(tmpE,aes(x=prolaris,y=oncotype))+
  geom_point(colour="tomato",size=3) +
  geom_point(data=tmpA,aes(x=prolaris,y=oncotype),colour="deepskyblue",size=3) +
  annotate("text",label="EAM = 0.54 (0.45 , 0.64)",x=119,y=-15,cex=5) +
  annotate("text",label="AAM = 0.63 (0.52 , 0.78)",x=119,y=-17,cex=5) +
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        axis.title = element_text(face="bold", size=15)) +
  geom_smooth(se=FALSE, colour="tomato", method="glm",fullrange=TRUE)+
  geom_smooth(data=tmpA,aes(x=prolaris,y=oncotype),se=FALSE, colour="deepskyblue", method="glm",fullrange=TRUE)+
  labs(x="Prolaris",y="Oncotype")

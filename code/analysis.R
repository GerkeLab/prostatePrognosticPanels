#####-----------------------------------------------------------------------------------------#####
# 01 - SETTING UP THE ENVIRONEMNT
#####-----------------------------------------------------------------------------------------#####

# import libraries used 
library(readr)
library(corrplot)
library(boot)
library(ggplot2)
library(reshape)
library(here)
library(ppcor)
library(viridis)

#####-----------------------------------------------------------------------------------------#####
# 02 - HELPER FUNCTION(S)
#####-----------------------------------------------------------------------------------------#####

# function used in calulcating confidence intervals for panel correlations 
## see section 07.02
corrFun <- function(data,indices){
  d <- data[indices,] # alows boot to select samples
  return(cor.test(d[,1],d[,2], method = "spearman")$estimate[[1]])
}

#####-----------------------------------------------------------------------------------------#####
# 03 - GENES IN EACH PANEL
#####-----------------------------------------------------------------------------------------#####

decipher_genes <- c("LASP1","IQGAP3","NFIB","S1PR4","THBS2","ANO7","PCDH7",
                    "MYBPC1","EPPK1","TSBP", "PBX1","NUSAP1","ZWILCH","UBE2C",
                    "CAMK2N1","RABGAP1","PCAT-32","PCAT-80","TNFRSF19","C6orf10")

prolaris_genes <- c("FOXM1","CDC20","CDKN3","CDC2","KIF11","KIAA0101","NUSAP1","CENPF","ASPM",
                    "BUB1B","RRM2","DLGAP5","BIRC5","KIF20A","PLK1","TOP2A","TK1","PBK","ASF1B",
                    "C18orf24","RAD54L","PTTG1","CDCA3","MCM10","PRC1","DTL","CEP55","RAD51",
                    "CENPM","CDCA8","ORC6L","SKA1","ORC6","CDK1") 

oncotypedx_genes <- c("AZGP1","KLK2","SRD5A2","FAM13C","FLNC","GSN","TPM2","GSTM2","TPX2","BGN",
                      "COL1A1","SFRP4")

#####-----------------------------------------------------------------------------------------#####
# 04 - DATA IMPORT 
#####-----------------------------------------------------------------------------------------#####

# single file containing clinical + expression data 
dat <- read_delim(here("data/nano.txt"), "\t", escape_double = FALSE, trim_ws = TRUE)
dat <- as.data.frame(dat)

dat$CAPRASn <- ifelse(dat$CAPRASgroup=='LOW',0,
                      ifelse(dat$CAPRASgroup=="HIGH",2,1))

#####-----------------------------------------------------------------------------------------#####
# 05 - DIFFERNETIAL XPRESSION - EAM VS AAM
#####-----------------------------------------------------------------------------------------#####

# calculate p value for mann whitney between races for each gene
raceMannAll <- apply(dat[,colnames(dat) %in% c(prolaris_genes,oncotypedx_genes,decipher_genes)],
                     2,function(x) wilcox.test(x~dat$Race)$p.value) 
raceMannAll <- as.data.frame(cbind(gene=colnames(dat[,colnames(dat) %in% c(prolaris_genes,oncotypedx_genes,decipher_genes)]),
                                   raceMannAll))
raceMannAll$raceMannAll <- as.numeric(as.character(raceMannAll$raceMannAll))
raceMannAll$gene <- as.character(raceMannAll$gene)

# annoate labels
raceMannAll$label <- ifelse(raceMannAll$raceMannAll<=0.001,"***",
                            ifelse(raceMannAll$raceMannAll<=0.01,"**",
                                   ifelse(raceMannAll$raceMannAll<=0.05,"*",NA)))

# calculate fold change
raceMannAll$fold <- NA

for(i in 1:length(raceMannAll$gene)){
  raceMannAll$fold[i] <- round((median(dat[dat$Race=="EAM",as.character(raceMannAll$gene[i])])/median(dat[dat$Race=="AAM",as.character(raceMannAll$gene[i])])),2)
}

# raceMannAll$fold <- round(log(raceMannAll$fold,2),2)
raceMannAll <- raceMannAll[!is.na(raceMannAll$label),]

melted_data <- dat[,colnames(dat) %in% c("ID","Race",prolaris_genes,oncotypedx_genes,decipher_genes)]
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
  scale_color_viridis_d()

### Figure 1B
ggplot(melted_data[melted_data$variable %in% oncotypedx_genes,],
       aes(x=variable, y=value, fill=as.factor(Race)))+
  scale_fill_manual(values=rep("white",length(oncotypedx_genes)))+
  geom_jitter(position=position_jitterdodge(jitter.width=.15, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.25)+
  geom_boxplot(outlier.shape=NA, alpha=0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                   face = c(rep("italic",4),rep("bold",3),
                                            rep("italic",3), "bold", "italic"),
                                   color = "black", size=10),
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
  scale_color_viridis_d(labels=c("AAM","EAM"))

### Figure 1C
ggplot(melted_data[melted_data$variable %in% decipher_genes,],
       aes(x=variable, y=value, fill=as.factor(Race)))+
  scale_fill_manual(values=rep("white",length(decipher_genes)))+
  geom_jitter(position=position_jitterdodge(jitter.width=.15, jitter.height=0, dodge.width=.75),
              aes(fill=as.factor(Race), col=as.factor(Race)),alpha=0.25)+
  geom_boxplot(outlier.shape=NA, alpha=0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   face = c(rep("bold",2), rep("italic",4),
                                            "bold", rep("italic",5),
                                            rep("bold",2), rep("italic",2),
                                            "bold", rep("italic",2)),
                                   color = "black", size=10),
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
  scale_color_viridis_d()

#####-----------------------------------------------------------------------------------------#####
# 06 - INTER-GENE CORRELATION (SEPERATELY IN EAM AND AAM) AND ICC
#####-----------------------------------------------------------------------------------------#####

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
         # cl.lim = c(-0.5,1), method="color", col= colorRampPalette(c("blue","white", "red"))(20),
         cl.lim = c(-0.5,1), method="color", col= viridis(20),
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
         cl.lim = c(-0.5,1), method="color", col= viridis(20),
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

#####-----------------------------------------------------------------------------------------#####
# Oncotype DX

## Calculate spearmans correlations rounded to one decimal place
corrAAM <- round(cor(dat[dat$Race=="AAM",colnames(dat) %in% oncotypedx_genes],method="spearman"),1)
corrEAM <- round(cor(dat[dat$Race=="EAM",colnames(dat) %in% oncotypedx_genes],method="spearman"),1)

### Figure 2B1
corrplot(corrEAM[c("TPX2","SFRP4","BGN","COL1A1","AZGP1","FAM13C","KLK2","SRD5A2","GSN","GSTM2","FLNC","TPM2"),
                 c("TPX2","SFRP4","BGN","COL1A1","AZGP1","FAM13C","KLK2","SRD5A2","GSN","GSTM2","FLNC","TPM2")],
         cl.lim = c(-0.5,1), method="color", col= viridis(20),
         tl.cex=1.5, tl.col = "black",cl.cex=0.6, cl.ratio = 0.1,addgrid.col = "lightgrey",
         family="sans", font=3)

### Figure 2B2
corrplot(corrAAM[c("TPX2","SFRP4","BGN","COL1A1","AZGP1","FAM13C","KLK2","SRD5A2","GSN","GSTM2","FLNC","TPM2"),
                 c("TPX2","SFRP4","BGN","COL1A1","AZGP1","FAM13C","KLK2","SRD5A2","GSN","GSTM2","FLNC","TPM2")],
         cl.lim = c(-0.5,1), method="color", col= viridis(20),
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

#####-----------------------------------------------------------------------------------------#####
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
         cl.lim = c(-0.5,1), method="color", col= viridis(20),
         tl.cex=1.3, tl.col = "black",cl.cex=0.6, cl.ratio = 0.1,addgrid.col = "lightgrey",
         family="sans",cl.pos="n",font=3)

### Figure 2C2
corrplot(corrAAM[c("PBX1","PCDH7","EPPK1","CAMK2N1","THBS2","ZWILCH","S1PR4","C6orf10","NUSAP1",
                   "UBE2C","IQGAP3","PCAT-32","PCAT-80","TNFRSF19","ANO7","LASP1","RABGAP1",
                   "MYBPC1","NFIB"),
                 c("PBX1","PCDH7","EPPK1","CAMK2N1","THBS2","ZWILCH","S1PR4","C6orf10","NUSAP1",
                   "UBE2C","IQGAP3","PCAT-32","PCAT-80","TNFRSF19","ANO7","LASP1","RABGAP1",
                   "MYBPC1","NFIB")],
         cl.lim = c(-0.5,1), method="color", col= viridis(20),
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

#####-----------------------------------------------------------------------------------------#####
# 07 - RISK SCORES - EAM VS AAM
#####-----------------------------------------------------------------------------------------#####
# Prolaris

## calculate risk score
prolarisGenes <- dat[,colnames(dat) %in% prolaris_genes]
gene_med <- apply(prolarisGenes,2,median)
prolarisGenes_centered <- prolarisGenes - gene_med
prolarisGenes_centered <- prolarisGenes_centered^2 # squaring the median centered expression values 

prolarisScore <- apply(prolarisGenes_centered,1,mean)
prolarisScore <- log2(prolarisScore)

risk_prolaris <- as.data.frame(cbind(id=dat$ID,
                                     score=prolarisScore,
                                     Race=dat$Race,
                                     capra=dat$CAPRASgroup,
                                     nccn=toupper(dat$nccn)))

risk_prolaris$score <- as.numeric(as.character(risk_prolaris$score))

wilcox.test(risk_prolaris$score~risk_prolaris$Race)

### Figure 3B 
B1 <- ggplot(risk_prolaris, aes(x=capra, y=score, fill=Race))+
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
  annotate("text",label=paste0("p = ",0.18),x="HIGH",y=0,cex=5, family="Calibri") + 
  scale_color_viridis_d()

### Figure 3B2 
B2 <- ggplot(risk_prolaris, aes(x=nccn, y=score, fill=Race))+
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
  labs(x="NCCN",y="Risk Score",title="Prolaris") +
  annotate("text",label=paste0("p = ",0.24),x="HIGH",y=0,cex=5, family="Calibri") + 
  scale_color_viridis_d()

#####-----------------------------------------------------------------------------------------#####
# oncotype dx

#  lower bound of expression values for TPX2 and SRDA5 were capped at 5 and 5.5 respectively
dat$TPX2_bounded <- ifelse(dat$TPX2 < 5, 5, dat$TPX2)
dat$SRD5A2_bounded <- ifelse(dat$SRD5A2 < 5.5, 5.5, dat$SRD5A2)

cellular_organization_module = dat$FLNC + dat$GSN + dat$TPM2 + dat$GSTM2
stromal_module = dat$BGN + dat$COL1A1 + dat$SFRP4
androgen_module = dat$FAM13C + dat$KLK2 + dat$SRD5A2_bounded + dat$AZGP1
proliferation_module = dat$TPX2_bounded

oncotypeScore = 0.735*stromal_module - 0.368*cellular_organization_module -
  0.352*androgen_module + 0.95*proliferation_module
risk_oncotype <- as.data.frame(cbind(id=dat$ID,
                                     score=oncotypeScore,
                                     Race=dat$Race,
                                     capra=dat$CAPRASgroup,
                                     nccn=toupper(dat$nccn)))
risk_oncotype$score <- as.numeric(as.character(risk_oncotype$score))

wilcox.test(risk_oncotype$score~risk_oncotype$Race)

## Figure 3C
C1 <- ggplot(risk_oncotype, aes(x=capra, y=score, fill=as.factor(Race)))+
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
  annotate("text",label=paste0("p = ",0.0009),x="HIGH",y=-10,cex=5, family="Calibri") + 
  scale_color_viridis_d()

## Figure 3C2
C2 <- ggplot(risk_oncotype, aes(x=nccn, y=score, fill=as.factor(Race)))+
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
  labs(x="NCCN",y="Risk Score",title="Oncotype DX") +
  annotate("text",label=paste0("p = ",0.01),x="HIGH",y=-10,cex=5, family="Calibri") + 
  scale_color_viridis_d()


#####-----------------------------------------------------------------------------------------#####
# Decipher

decipherGenes <- dat[,colnames(dat) %in% decipher_genes]
gene_med <- apply(decipherGenes,2,median)
decipherGenes_centered <- decipherGenes - gene_med

cluster1 <- hclust(dist(t(decipherGenes)),method="centroid")

over <- c("CAMK2N1","EPPK1","IQGAP3","LASP1","NFIB","NUSAP1","PBX1","S1PR4","THBS2","UBE2C",
          "ZWILCH")

under <- c("ANO7","C6orf10","MYBPC1","PCDH7","RABGAP1","TNFRSF19")

# average of the log2 normalized values for the 9 over-expressed targets
c1 <- apply(decipherGenes[,c("CAMK2N1","EPPK1","IQGAP3","LASP1","NFIB","NUSAP1","PBX1",
                             "S1PR4","THBS2","UBE2C","ZWILCH")],1,mean)

# average of the log2 normalized values for the 9 under-expressed targets
c2 <- apply(decipherGenes[,c("ANO7","C6orf10","MYBPC1","PCDH7","RABGAP1","TNFRSF19")],1,mean)

decipherScore <- c1-c2

risk_decipher <- as.data.frame(cbind(id=dat$ID,
                                     score=decipherScore,
                                     Race=dat$Race,
                                     capra=dat$CAPRASgroup,
                                     nccn=toupper(dat$nccn)))
risk_decipher$score <- as.numeric(as.character(risk_decipher$score))

wilcox.test(risk_decipher$score~risk_decipher$Race)

### Figure 3D
D1 <- ggplot(risk_decipher, aes(x=capra, y=score, fill=as.factor(Race)))+
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
  annotate("text",label=paste0("p = ",0.38),x="HIGH",y=-2.2,cex=5, family="Calibri") + 
  scale_color_viridis_d()

### Figure 3D2
D2 <- ggplot(risk_decipher, aes(x=nccn, y=score, fill=as.factor(Race)))+
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
  labs(x="NCCN",y="Risk Score",title="Decipher") +
  annotate("text",label=paste0("p = ",0.30),x="HIGH",y=-2.2,cex=5, family="Calibri") + 
  scale_color_viridis_d()

#####-----------------------------------------------------------------------------------------#####

### Figure 3A

tmp <- as.data.frame(cbind(Race=c("AAM","AAM","AAM","EAM","EAM","EAM"),
                           capra=c("LOW","INTERMEDIATE","HIGH","LOW","INTERMEDIATE","HIGH"),
                           perc=c(65,29,5,60,30,10)))

A1 <- ggplot(tmp, aes(fill=Race, x=capra)) +
  geom_col(aes(y=as.numeric(as.character(perc))),position = "dodge") + 
  scale_x_discrete(limits = c("LOW","INTERMEDIATE","HIGH")) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        # text = element_text(family="Calibri"),
        axis.title = element_text(face="bold", size=15)) +
  labs(x="CAPRA-S",y="%") +
  guides(fill=guide_legend(title="Race")) +
  scale_fill_viridis_d()

### Figure 3A2

tmp <- as.data.frame(cbind(Race=c("AAM","AAM","AAM","EAM","EAM","EAM"),
                           nccn=c("LOW","INTERMEDIATE","HIGH","LOW","INTERMEDIATE","HIGH"),
                           perc=c(4,77,19,9,54,37)))

A2 <- ggplot(tmp, aes(fill=Race, x=nccn)) +
  geom_col(aes(y=as.numeric(as.character(perc))),position = "dodge") + 
  scale_x_discrete(limits = c("LOW","INTERMEDIATE","HIGH")) + 
  theme(axis.text.x = element_text(color = "black", size=10),
        axis.text.y = element_text(color = "black", size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), 
        plot.title = element_text(hjust = 0.5, size=20, face="bold"),
        # text = element_text(family="Calibri"),
        axis.title = element_text(face="bold", size=15)) +
  labs(x="NCCN",y="%") +
  guides(fill=guide_legend(title="Race")) +
  scale_fill_viridis_d()

# pdf("/Volumes/Lab_Gerke/prostateWorkGroup/teamScienceGenes/Panelpaper/figures/riskScoresNCCN.pdf",
#     width = 8, height = 11)
gridExtra::grid.arrange(A1,A2,B1,B2,C1,C2,D1,D2, ncol=2)
# dev.off()

rm(A1,A2,B1,B2,C1,C2,D1,D2)

#####-----------------------------------------------------------------------------------------#####
# 07.02 - Panel risk score correlations 
#####-----------------------------------------------------------------------------------------#####

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
boot_results <- boot(data=tmpE[,1:2], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")
boot_results <- boot(data=tmpA[,1:2], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

# Figure 4C
ggplot(tmpE,aes(x=decipher,y=prolaris))+
  geom_point(colour="tomato",size=3) +
  geom_point(data=tmpA,aes(x=decipher,y=prolaris),colour="deepskyblue",size=3) +
  annotate("text",label="EAM = 0.30 (0.19 , 0.42)",x=0.8,y=0.2,cex=5) +
  annotate("text",label="AAM = 0.12 (-0.10 , 0.33)",x=0.8,y=0,cex=5) +
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

#####-----------------------------------------------------------------------------------------#####

# decipher vs oncotype
cor.test(tmpE$decipher,tmpE$oncotype, method = "spearman")
cor.test(tmpA$decipher,tmpA$oncotype, method = "spearman")

# bootstrap for CI
boot_results <- boot(data=tmpE[,c(1,3)], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")
boot_results <- boot(data=tmpA[,c(1,3)], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

# Fifure 4B
ggplot(tmpE,aes(x=decipher,y=oncotype))+
  geom_point(colour="tomato",size=3) +
  geom_point(data=tmpA,aes(x=decipher,y=oncotype),colour="deepskyblue",size=3) +
  annotate("text",label="EAM = 0.66 (0.59 , 0.75)",x=0.8,y=-9,cex=5) +
  annotate("text",label="AAM = 0.61 (0.50 , 0.76)",x=0.8,y=-10,cex=5) +
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


#####-----------------------------------------------------------------------------------------#####

# prolaris vs oncotype
cor.test(tmpE$prolaris,tmpE$oncotype, method = "spearman")
cor.test(tmpA$prolaris,tmpA$oncotype, method = "spearman")

boot_results <- boot(data=tmpE[,2:3], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")
boot_results <- boot(data=tmpA[,2:3], statistic=corrFun,R=1000)
boot.ci(boot_results, type="basic")

ggplot(tmpE,aes(x=prolaris,y=oncotype))+
  geom_point(colour="tomato",size=3) +
  geom_point(data=tmpA,aes(x=prolaris,y=oncotype),colour="deepskyblue",size=3) +
  annotate("text",label="EAM = 0.26 (0.14 , 0.38)",x=0.2,y=2,cex=5) +
  annotate("text",label="AAM = 0.17 (-0.02 , 0.37)",x=0.2,y=1,cex=5) +
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

#####-----------------------------------------------------------------------------------------#####
# 07.03 - Risk scores and clinical features
#####-----------------------------------------------------------------------------------------#####

risk_decipher <- data.frame(id=dat$ID,
                            score=decipherScore,
                            Race=dat$Race,
                            capra=dat$CAPRASn,
                            nccn=dat$nccn)
risk_decipher <- risk_decipher %>%
  mutate(score = as.numeric(as.character(score))) %>%
  mutate(nccn = as.numeric(factor(nccn, levels = rev(levels(nccn)))))

cor.test(risk_decipher$score, risk_decipher$capra)
cor.test(risk_decipher[risk_decipher$Race=="EAM",]$score, risk_decipher[risk_decipher$Race=="EAM",]$capra) 
cor.test(risk_decipher[risk_decipher$Race=="AAM",]$score, risk_decipher[risk_decipher$Race=="AAM",]$capra) 

cor.test(risk_decipher$score, risk_decipher$nccn)
cor.test(risk_decipher[risk_decipher$Race=="EAM",]$score, risk_decipher[risk_decipher$Race=="EAM",]$nccn) 
cor.test(risk_decipher[risk_decipher$Race=="AAM",]$score, risk_decipher[risk_decipher$Race=="AAM",]$nccn) 

risk_prolaris <- data.frame(id=dat$ID,
                            score=prolarisScore,
                            Race=dat$Race,
                            capra=dat$CAPRASn,
                            nccn=dat$nccn)
risk_prolaris <- risk_prolaris %>%
  mutate(score = as.numeric(as.character(score))) %>%
  mutate(nccn = as.numeric(factor(nccn, levels = rev(levels(nccn)))))

cor.test(risk_prolaris$score, risk_prolaris$capra)
cor.test(risk_prolaris[risk_prolaris$Race=="EAM",]$score, risk_prolaris[risk_prolaris$Race=="EAM",]$capra) 
cor.test(risk_prolaris[risk_prolaris$Race=="AAM",]$score, risk_prolaris[risk_prolaris$Race=="AAM",]$capra) 

cor.test(risk_prolaris$score, risk_prolaris$nccn)
cor.test(risk_prolaris[risk_prolaris$Race=="EAM",]$score, risk_prolaris[risk_prolaris$Race=="EAM",]$nccn) 
cor.test(risk_prolaris[risk_prolaris$Race=="AAM",]$score, risk_prolaris[risk_prolaris$Race=="AAM",]$nccn) 

risk_oncotype <- data.frame(id=dat$ID,
                            score=oncotypeScore,
                            Race=dat$Race,
                            capra=dat$CAPRASn,
                            nccn=dat$nccn)
risk_oncotype <- risk_oncotype %>%
  mutate(score = as.numeric(as.character(score))) %>%
  mutate(nccn = as.numeric(factor(nccn, levels = rev(levels(nccn)))))

cor.test(risk_oncotype$score, risk_oncotype$capra)
cor.test(risk_oncotype[risk_oncotype$Race=="EAM",]$score, risk_oncotype[risk_oncotype$Race=="EAM",]$capra) 
cor.test(risk_oncotype[risk_oncotype$Race=="AAM",]$score, risk_oncotype[risk_oncotype$Race=="AAM",]$capra) 

cor.test(risk_oncotype$score, risk_oncotype$nccn)
cor.test(risk_oncotype[risk_oncotype$Race=="EAM",]$score, risk_oncotype[risk_oncotype$Race=="EAM",]$nccn) 
cor.test(risk_oncotype[risk_oncotype$Race=="AAM",]$score, risk_oncotype[risk_oncotype$Race=="AAM",]$nccn) 

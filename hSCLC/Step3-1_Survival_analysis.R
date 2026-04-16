library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(survival)
library(survminer)
library(readxl)
library(readxlsb)

# Read log2TPM+1 data from the Supplementary Table S1D of Peng Zhang et al. Cell 2024
log2TPM <- readxlsb::read_xlsb("Table S1.xlsb",sheet = "Table S1D",col_names = T) %>% as.data.frame()
rownames(log2TPM) <- log2TPM$Gene
log2TPM <- log2TPM[,-1] %>% as.matrix()
log2TPM[is.na(log2TPM)] <- 0
log2TPM.scaled <- scale(t(log2TPM),center = T,scale = T) %>% t()

# Read Survaival info from Supplementary Table S1A of Peng Zhang et al. Cell 2024
surv <- readxlsb::read_xlsb("Table S1.xlsb",sheet = "Table S1A",col_names = T) %>% as.data.frame()
surv <- surv[,c("Sample.ID","Histologic.type","Status.","Survial..months.")]
surv$Survial..months. <- as.numeric(surv$Survial..months.)
surv$Sample.ID <- paste0("T",surv$Sample.ID)
surv$OS_status <- ifelse(surv$Status.=="alive",0,1)

# Read the N-E score mouse reference data, supplied in Supplementary Table
NE_marker50 <- readxl::read_xlsx("NE_markers_ref_top50_top25.xlsx",sheet = 1) %>% as.data.frame()
NE_marker50 <- NE_marker50[which(NE_marker50$humanGene %in% rownames(log2TPM)),]

# Calculate the N-E score of patients
NE_score1 <- c()
pb <- txtProgressBar(min = 1,max = ncol(log2TPM),style = 3)
for (i in 1:ncol(log2TPM)) {
  tmp <- log2TPM[NE_marker50$humanGene,i]
  cor_N <- cor.test(tmp,NE_marker50$N_expression,method = "spearman")$estimate
  cor_E <- cor.test(tmp,NE_marker50$E_expression,method = "spearman")$estimate
  NE_score1[i] <- (cor_N-cor_E)/2
  
  setTxtProgressBar(pb,i)
  rm(tmp,cor_N,cor_E)
}
NE_score <- data.frame(sample=colnames(log2TPM),NE_spearman50=NE_score1)
plotdata <- merge(surv,NE_score,by.x="Sample.ID",by.y="sample")
plotdata <- plotdata[which(plotdata$Histologic.type=="small cell lung cancer"),]

# Define N-E score group
cutoff <- surv_cutpoint(plotdata,time = "Survial..months.",event = "OS_status",variables = "NE_spearman50",minprop = 0.1)
plotdata$NE50_group <- ifelse(plotdata$NE_spearman50>=cutoff$cutpoint[1,1],"high","low")

# Survival analysis
pdf("Result/Survival_all.pdf",height = 5,width = 6)
fit <- survfit(Surv(time = `Survial..months.`,event = OS_status) ~ NE50_group, data = plotdata)
ggsurvplot(fit,conf.int = F,risk.table = T,risk.table.col = "strata",pval = T,
           break.time = 10,surv.median.line = "hv",title="All_by_NE_spearman_score")
dev.off()



options(stringsAsFactors = FALSE)

source("JBC_DNM_func.r")

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggExtra)
library(rms)

# =========================================================
# Figure 2B. Mutation-type-specific parental age effects
# =========================================================

covar_list = c("age_mat","age_pat")
covar_class = rep("c",2)
covar_name = c("Maternal age","Paternal age")

pheno_Ps = paste0("dSNVP_",c("TG","TC","TA","CA","CG","CT_CpG","CT_nCpG"))
pheno_Ms = paste0("dSNVM_",c("TG","TC","TA","CA","CG","CT_CpG","CT_nCpG"))

Ps_result <- rbindlist(lapply(pheno_Ps, function(p) {
    res_art <- generate_multi_linear(dnm_jbc_input,covar_list,covar_class,covar_name,p)
    res_art$pheno <- p
	res_art$mutType <- sub("dSNVP_", "", p)
	res_art$mutType[res_art$mutType == "CT_CpG"]  <- "C>T (CpG sites)"
	res_art$mutType[res_art$mutType == "CT_nCpG"] <- "C>T (Other)"
	res_art$mutType <- gsub("^([ACGT])([ACGT])$", "\\1>\\2", res_art$mutType)
    res_art[2,] 
}),fill = TRUE)

Ms_result <- rbindlist(lapply(pheno_Ms, function(p) {
    res_art <- generate_multi_linear(dnm_jbc_input,covar_list,covar_class,covar_name,p)
    res_art$pheno <- p
	res_art$mutType <- sub("dSNVM_", "", p)
	res_art$mutType[res_art$mutType == "CT_CpG"]  <- "C>T (CpG sites)"
	res_art$mutType[res_art$mutType == "CT_nCpG"] <- "C>T (Other)"
	res_art$mutType <- gsub("^([ACGT])([ACGT])$", "\\1>\\2", res_art$mutType)
    res_art[1,] 
}),fill = TRUE)


figure2b_data <- cbind(Ps_result[,c(1,5:6,10)],Ms_result[,c(1,5:6)])
colnames(figure2b_data)[1:3] <- paste0(colnames(figure2b_data)[1:3],"_P")
colnames(figure2b_data)[5:7] <- paste0(colnames(figure2b_data)[5:7],"_M")

figure2b_data[,c(2:3,6:7)] <- lapply(figure2b_data[,c(2:3,6:7)], as.numeric)

ggplot(figure2b_data, aes(Beta_M,Beta_P,col=mutType)) +
	theme_classic() +
	geom_point(position = position_dodge(0.4), size=2)+
	geom_errorbar(aes(ymin=Beta_P-1.96*SE_P, ymax=Beta_P+1.96*SE_P),width=0) +
	geom_errorbarh(aes(xmin=Beta_M-1.96*SE_M, xmax=Beta_M+1.96*SE_M), height=0) +
	labs(x="Maternal age effect", y="Paternal age effect") +
	xlim(-0.05,0.15) +
	scale_colour_manual(
		values=c("#e1061f","#FBB255","#1977b7","#38a750","#8e498f","#ef801a","#ea83b0"),
		breaks=c("C>A","C>G","C>T (CpG sites)","C>T (Other)","T>A","T>C","T>G"),
		labels=c("C>A","C>G","C>T (CpG sites)","C>T (Other)","T>A","T>C","T>G")
	) +
	geom_smooth(method="lm", se=FALSE, colour="black")  +
	theme(
		axis.text.y=element_text(size=10,color="#030303",face="bold"),
		axis.title.y=element_text(size=10,face="bold"),
		axis.text.x=element_text(face="bold",color="black",size=10),
		axis.title.x=element_text(face="bold",color="black",size=10),
		panel.grid.major=element_line(colour=NA),
		panel.grid.minor=element_line(colour=NA)
	)

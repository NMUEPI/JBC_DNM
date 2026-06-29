library(MatchIt)
library(optmatch)
library(ggplot2)
library(dplyr)
library(rms)

DNM_mut_all_final <- read.delim("dnm_final_list_forcheck.tsv",h=T)
dnm_final_foranalysis <- read.delim("dnm_final_data_forcheck.tsv",h=T)
dnm_submit_data<-read.delim("dnm_jbc_input_20260625.txt",h=T)

dnm_final_foranalysis_singleton <- subset(dnm_final_foranalysis,fetus_n==1)
dnm_final_foranalysis_ART <- subset(dnm_final_foranalysis,group_add=="ART")
dnm_final_foranalysis_SP <- subset(dnm_final_foranalysis,group_add=="SP")

# =========================================================
# Part 1 General description of the study population and de novo mutations
# =========================================================

################ Supplementary Table 1
table1 <- table_baseline_final(subset(dnm_final_foranalysis,!duplicated(cohortid)),"group_add")

nrow(RF_info_foranalysis)/8328/mean(dnm_final_foranalysis$callableregion_length)/2

################ Figure 1B

DNM_list <- list(
    FPOG      = subset(DNM_mut_all_final, Source2 %in% c("paternal","maternal")),
    EPZM      = subset(DNM_mut_all_final, Source2 == "mosaic"),
    Paternal  = subset(DNM_mut_all_final, Source2 == "paternal"),
    Maternal  = subset(DNM_mut_all_final, Source2 == "maternal")
)

mutation_type <- c("C>A","C>G","C>T(CpG)","C>T(Other)","T>G","T>C","T>A")

figure1b_result <- data.frame(type = mutation_type)

for(i in seq_along(mutation_type)){

    mt <- mutation_type[i]

    ## FPOG vs EPZM

    tab1 <- matrix(c(
        sum(DNM_list$FPOG$nchanget==mt),
        sum(DNM_list$FPOG$nchanget!=mt),
        sum(DNM_list$EPZM$nchanget==mt),
        sum(DNM_list$EPZM$nchanget!=mt)
    ),2)

    test1 <- fisher.test(tab1)

    ## Paternal vs Maternal

    tab2 <- matrix(c(
        sum(DNM_list$Paternal$nchanget==mt),
        sum(DNM_list$Paternal$nchanget!=mt),
        sum(DNM_list$Maternal$nchanget==mt),
        sum(DNM_list$Maternal$nchanget!=mt)
    ),2)

    test2 <- fisher.test(tab2)

    figure1b_result$P1[i] <- test1$p.value
    figure1b_result$P2[i] <- test2$p.value

    figure1b_result$prop_P[i] <- mean(DNM_list$Paternal$nchanget==mt)
    figure1b_result$prop_M[i] <- mean(DNM_list$Maternal$nchanget==mt)
    figure1b_result$prop_EPZM[i] <- mean(DNM_list$EPZM$nchanget==mt)

}

figure1b_result$padj1 <- p.adjust(figure1b_result$P1,"fdr")
figure1b_result$padj2 <- p.adjust(figure1b_result$P2,"fdr")

figure1b_data <- bind_rows(
    data.frame(
        muttype = figure1b_result$type,
        prop = figure1b_result$prop_P,
        Source = "Paternal dnSNV"),
    data.frame(
        muttype = figure1b_result$type,
        prop = figure1b_result$prop_M,
        Source = "Maternal dnSNV"),
    data.frame(
        muttype = figure1b_result$type,
        prop = figure1b_result$prop_EPZM,
        Source = "EPZM")
)

pdf("Figure1B.pdf",width=10,height=4)
figure1b_data$muttype <- factor(figure1b_data$muttype,levels=c("C>A","C>G","C>T(CpG)","C>T(Other)","T>A","T>C","T>G"))
figure1b_data$Source <- factor(figure1b_data$Source,levels=c("Maternal dnSNV","Paternal dnSNV","EPZM"))

f <- ggplot(figure1b_data, aes(muttype, prop,fill=Source))+geom_bar(stat = "identity", position=position_dodge(),color="grey",width=0.7)+theme_bw()+theme_classic()+theme(axis.text.y=element_text(size=12,color="#030303",face="bold"),axis.title.y=element_text(size=12,face="bold"),axis.text.x=element_text(face="bold",color="black",size=12),axis.title.x=element_text(face="bold",color="black",size=12),panel.grid.major=element_line(colour=NA))+ theme(legend.title = element_blank())+xlab(NULL)+ylab("Mutation Proportion")+ scale_fill_manual(values=c("#CC1B21","#0066BD","#FBB255"),breaks=c("Maternal dnSNV","Paternal dnSNV","EPZM"),labels=c("Maternal dnSNV","Paternal dnSNV","EPZM"))

f + annotate("rect", xmin = 0.8, xmax = 1.0, ymin = 0.13, ymax =0.13, alpha=1,colour = "black")+
    annotate("rect", xmin = 0.8, xmax = 0.8, ymin = 0.11, ymax =0.13, alpha=1, colour = "black")+
    annotate("rect", xmin = 1.0, xmax = 1.0, ymin = 0.12, ymax =0.13, alpha=1, colour = "black")+

    annotate("rect", xmin = 0.9, xmax = 1.2, ymin = 0.15, ymax =0.15, alpha=1,colour = "black") +
    annotate("rect", xmin = 0.9, xmax = 0.9, ymin = 0.135, ymax =0.15, alpha=1, colour = "black")+
    annotate("rect", xmin = 1.2, xmax = 1.2, ymin = 0.13, ymax =0.15, alpha=1, colour = "black")+

    annotate("rect", xmin = 1.8, xmax = 2.0, ymin = 0.12, ymax =0.12, alpha=1,colour = "black") +
    annotate("rect", xmin = 1.8, xmax = 1.8, ymin = 0.10, ymax =0.12, alpha=1,colour = "black") +
    annotate("rect", xmin = 2.0, xmax = 2.0, ymin = 0.11, ymax =0.12, alpha=1,colour = "black") +

    annotate("rect", xmin = 2.8, xmax = 3.0, ymin = 0.11, ymax =0.11, alpha=1,colour = "black") +
    annotate("rect", xmin = 2.8, xmax = 2.8, ymin = 0.09, ymax =0.11, alpha=1,colour = "black") +
    annotate("rect", xmin = 3.0, xmax = 3.0, ymin = 0.10, ymax =0.11, alpha=1,colour = "black") +

    annotate("rect", xmin = 3.8, xmax = 4.0, ymin = 0.33, ymax =0.33, alpha=1,colour = "black") +
    annotate("rect", xmin = 3.8, xmax = 3.8, ymin = 0.32, ymax =0.33, alpha=1,colour = "black") +
    annotate("rect", xmin = 4.0, xmax = 4.0, ymin = 0.265, ymax =0.33, alpha=1,colour = "black") +

    annotate("rect", xmin = 3.9, xmax = 4.2, ymin = 0.35, ymax =0.35, alpha=1,colour = "black") +
    annotate("rect", xmin = 3.9, xmax = 3.9, ymin = 0.34, ymax =0.35, alpha=1,colour = "black") +
    annotate("rect", xmin = 4.2, xmax = 4.2, ymin = 0.30, ymax =0.35, alpha=1,colour = "black") +
    
    annotate("rect", xmin = 4.8, xmax = 5.0, ymin = 0.08, ymax =0.08, alpha=1,colour = "black") +
    annotate("rect", xmin = 4.9, xmax = 5.2, ymin = 0.10, ymax =0.10, alpha=1,colour = "black") +
    annotate("rect", xmin = 4.9, xmax = 4.9, ymin = 0.08, ymax =0.10, alpha=1,colour = "black") +
    annotate("rect", xmin = 5.2, xmax = 5.2, ymin = 0.09, ymax =0.10, alpha=1,colour = "black") +

    annotate("rect", xmin = 5.8, xmax = 6.0, ymin = 0.31, ymax =0.31, alpha=1,colour = "black") +
    annotate("rect", xmin = 5.8, xmax = 5.8, ymin = 0.29, ymax =0.31, alpha=1,colour = "black") +
    annotate("rect", xmin = 6.0, xmax = 6.0, ymin = 0.30, ymax =0.31, alpha=1,colour = "black") +

    annotate("rect", xmin = 5.9, xmax = 6.2, ymin = 0.32, ymax =0.32, alpha=1,colour = "black") +
    annotate("rect", xmin = 5.9, xmax = 5.9, ymin = 0.31, ymax =0.32, alpha=1,colour = "black") +
    annotate("rect", xmin = 6.2, xmax = 6.2, ymin = 0.255, ymax =0.32, alpha=1,colour = "black") +

    annotate("rect", xmin = 6.8, xmax = 7.0, ymin = 0.09, ymax =0.09, alpha=1,colour = "black") +
    annotate("rect", xmin = 6.9, xmax = 7.2, ymin = 0.105, ymax =0.105, alpha=1,colour = "black") +
    annotate("rect", xmin = 6.9, xmax = 6.9, ymin = 0.09, ymax =0.105, alpha=1,colour = "black") +
    annotate("rect", xmin = 7.2, xmax = 7.2, ymin = 0.09, ymax =0.105, alpha=1,colour = "black")
dev.off()

################ Figure 1C

tricontext_type <- unique(DNM_mut_all_final$tricontext)

figure1c_result <- data.frame(muttype = tricontext_type)

for(i in seq_along(tricontext_type)){

    mt <- tricontext_type[i]

    ## FPOG vs EPZM

    tab1 <- matrix(c(
        sum(DNM_list$EPZM$tricontext==mt),
        sum(DNM_list$EPZM$tricontext!=mt),
        sum(DNM_list$FPOG$tricontext==mt),
        sum(DNM_list$FPOG$tricontext!=mt)
    ),2)

    test1 <- fisher.test(tab1)

    ## Paternal vs Maternal

    tab2 <- matrix(c(
        sum(DNM_list$Paternal$tricontext==mt),
        sum(DNM_list$Paternal$tricontext!=mt),
        sum(DNM_list$Maternal$tricontext==mt),
        sum(DNM_list$Maternal$tricontext!=mt)
    ),2)

    test2 <- fisher.test(tab2)

	figure1c_result$ER1[i] <- test1$estimate
    figure1c_result$ER2[i] <- test2$estimate
	
    figure1c_result$P1[i] <- test1$p.value
    figure1c_result$P2[i] <- test2$p.value

    figure1c_result$prop_P[i] <- mean(DNM_list$Paternal$tricontext==mt)
    figure1c_result$prop_M[i] <- mean(DNM_list$Maternal$tricontext==mt)
    figure1c_result$prop_EPZM[i] <- mean(DNM_list$EPZM$tricontext==mt)

}

mut_anno <- unique(DNM_mut_all_final[,c("tricontext","std.mutcat")])

figure1c_d1 <- rbind(data.frame(
        muttype = figure1c_result$muttype,
        ER      = figure1c_result$ER1,
        P       = figure1c_result$P1,
        Source  = "POG vs EPZM"),
    data.frame(
        muttype = figure1c_result$muttype,
        ER      = figure1c_result$ER2,
        P       = figure1c_result$P2,
        Source  = "Paternal vs Maternal")
)

figure1c_d1 <- merge(figure1c_d1,mut_anno, by.x="muttype", by.y="tricontext")

figure1c_d1$base1 <- sub("\\[.*","",figure1c_d1$muttype)
figure1c_d1$base2 <- sub(".*\\]","",figure1c_d1$muttype)

figure1c_d1$basenear <- factor(
    paste0(figure1c_d1$base1,"-",figure1c_d1$base2),
    levels=c("A-A","A-C","A-G","A-T",
             "C-A","C-C","C-G","C-T",
             "G-A","G-C","G-G","G-T",
             "T-A","T-C","T-G","T-T")
)

figure1c_d1$std.mutcat <- factor(
    figure1c_d1$std.mutcat,
    levels=c("C>A","C>G","C>T","T>A","T>C","T>G")
)

figure1c_d1$anno <- ifelse(figure1c_d1$P<0.05/96/2,"***","")

figure1c_d2 <- rbind(
    data.frame(
        muttype = figure1c_result$muttype,
        prop = figure1c_result$prop_P,
        Source = "Paternal dnSNV"),
    data.frame(
        muttype = figure1c_result$muttype,
        prop = figure1c_result$prop_M,
        Source = "Maternal dnSNV"),
    data.frame(
        muttype = figure1c_result$muttype,
        prop = figure1c_result$prop_EPZM,
        Source = "EPZM")
)

figure1c_d2 <- merge(
    figure1c_d2,
    mut_anno,
    by.x="muttype",
    by.y="tricontext"
)

figure1c_d2$base1 <- sub("\\[.*","",figure1c_d2$muttype)
figure1c_d2$base2 <- sub(".*\\]","",figure1c_d2$muttype)

figure1c_d2$basenear <- factor(
    paste0(figure1c_d2$base1,"-",figure1c_d2$base2),
    levels=c("A-A","A-C","A-G","A-T",
             "C-A","C-C","C-G","C-T",
             "G-A","G-C","G-G","G-T",
             "T-A","T-C","T-G","T-T")
)

figure1c_d2$std.mutcat <- factor(
    figure1c_d2$std.mutcat,
    levels=c("C>A","C>G","C>T","T>A","T>C","T>G")
)

library(cowplot)

color = c("#FEE07F","#E42927","#ECC7C5","#A2CF63","#507CB7","#91C0DB")
figure1c_d2$Source <- factor(figure1c_d2$Source,levels=c("Maternal dnSNV","Paternal dnSNV","EPZM"))

pdf("Figure1C.pdf",width = 13,height=4)
g1_up <- ggplot(figure1c_d2, aes(basenear, prop,fill=std.mutcat))+
		geom_bar(stat = "identity", position=position_dodge(),width=0.7)+
		theme_bw()+theme_classic()+
		theme(axis.text.y=element_text(size=10,color="#030303",face="bold"),axis.title.y=element_text(size=12,face="bold"),axis.text.x=element_text(color="black",size=10,angle=90, hjust = 1),axis.title.x=element_text(face="bold",color="black",size=12),panel.grid.major=element_line(colour=NA))+ 
		theme(legend.position = "none")+
		theme(strip.text.x = element_text(face = "bold", color = "white", hjust = 0, size = 15), strip.background.x = element_rect(fill = color, linetype = "solid",color = "white"), strip.background.y = element_rect(color = "white"))+
		xlab(NULL)+ylab("Mutation Proportion")+ 
		scale_fill_manual(values = color,breaks=c("C>A","C>G","C>T","T>A","T>C","T>G"),labels=c("C>A","C>G","C>T","T>A","T>C","T>G"))+facet_grid(Source~std.mutcat)
g1_down <- ggplot(figure1c_d1, aes(x=basenear, y=Source, fill=ER)) + geom_tile(color="white") + 
		theme_bw() + 
		theme(axis.text.y=element_text(color="black",size=10,hjust = 1),axis.title.y=element_text(size=12,face="bold"),axis.text.x=element_text(color="black",size=10,angle=90, hjust = 1),axis.title.x=element_text(face="bold",color="black",size=12),panel.grid.major=element_line(colour=NA)) +
		xlab("") + ylab("Mutation Proportion") +
		geom_text(label=figure1c_d1$anno, size=2) + 
		scale_fill_distiller(palette = "Spectral")+facet_grid(.~std.mutcat)

plot_grid(g1_up, g1_down,ncol=1, align='v', axis="lr",rel_heights=c(2,1))
dev.off()


########### Figure 1D

####
library(ggbreak)

num_func_epzm <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Functional" & Source2=="mosaic"))
num_func_pog <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Functional" & Source2=="paternal"))
num_func_mog <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Functional" & Source2=="maternal"))
num_func_pmog <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Functional" & Source2%in%c("maternal","paternal")))

num_NG_epzm <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Gene near" & Source2=="mosaic"))
num_NG_pog <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Gene near" & Source2=="paternal"))
num_NG_mog <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Gene near" & Source2=="maternal"))
num_NG_pmog <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Gene near" & Source2%in%c("maternal","paternal")))

num_IG_epzm <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Intergenic" & Source2=="mosaic"))
num_IG_pog <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Intergenic" & Source2=="paternal"))
num_IG_mog <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Intergenic" & Source2=="maternal"))
num_IG_pmog <- nrow(subset(DNM_mut_all_final,impact_gpCR=="Intergenic" & Source2%in%c("maternal","paternal")))

num_epzm <- nrow(subset(DNM_mut_all_final,Source2=="mosaic"))
num_pog <- nrow(subset(DNM_mut_all_final,Source2=="paternal"))
num_mog <- nrow(subset(DNM_mut_all_final,Source2=="maternal"))
num_pmog <- nrow(subset(DNM_mut_all_final,Source2%in%c("maternal","paternal")))

figure4_data <- rbind(	data.frame(muttype="Functional",prop=num_func_pog/num_pog*100,Source="Paternal DNM"),
						data.frame(muttype="Functional",prop=num_func_mog/num_mog*100,Source="Maternal DNM"),
						 data.frame(muttype="Functional",prop=num_func_epzm/num_epzm*100,Source="EPZM"),
						data.frame(muttype="Gene proximal",prop=num_NG_pog/num_pog*100,Source="Paternal DNM"),
						data.frame(muttype="Gene proximal",prop=num_NG_mog/num_mog*100,Source="Maternal DNM"),
						 data.frame(muttype="Gene proximal",prop=num_NG_epzm/num_epzm*100,Source="EPZM"),
						data.frame(muttype="Intergenic",prop=num_IG_pog/num_pog*100,Source="Paternal DNM"),
						data.frame(muttype="Intergenic",prop=num_IG_mog/num_mog*100,Source="Maternal DNM"),
						 data.frame(muttype="Intergenic",prop=num_IG_epzm/num_epzm*100,Source="EPZM"))

figure4_data$Source <- factor(figure4_data$Source,levels=c("Maternal DNM","Paternal DNM","EPZM"))
figure4_data$muttype <- factor(figure4_data$muttype,levels=c("Functional","Gene proximal","Intergenic"))

pdf("Figure1D.pdf",width = 8,height=4)
ggplot(figure4_data, aes(muttype, prop,fill=Source))+g
		eom_bar(stat = "identity", position=position_dodge(),color="grey",width=0.7)+
		theme_bw()+
		theme_classic()+
		theme(axis.text.y=element_text(size=12,color="#030303",face="bold"),
			axis.title.y=element_text(size=12,face="bold"),
			axis.text.x=element_text(face="bold",color="black",size=12),
			axis.title.x=element_text(face="bold",color="black",size=12),
			panel.grid.major=element_line(colour=NA))+ 
		theme(legend.title = element_blank())+ 
		xlab(NULL)+ylab("Mutation Proportion (%)")+ 
		scale_fill_manual(values=c("#CC1B21","#0066BD","#FBB255"),
			breaks=c("Maternal DNM","Paternal DNM","EPZM"),
			labels=c("Maternal DNM","Paternal DNM","EPZM"))
dev.off()

################################################################################################################
######### Part 2: age-related

## Figure 2B

figure2_alltype <- rbindlist(mclapply(1:length(phenoname_all),function(i) {
	print(i)
	figure2 <- univariable_model_onlyage(dnm_final_foranalysis,phenoname_all[i],"glmm")
	figure2$MutType <- phenoname_all[i]
	return(data.frame(figure2))
},mc.cores=40))
write.csv(figure2_alltype,paste0(path,"/tableS4.csv"),row.names=F)

figure2_alltype <- rbindlist(mclapply(1:length(phenoname_all),function(i) {
	print(i)
	figure2 <- univariable_model_onlyage(dnm_final_foranalysis_singleton,phenoname_all[i],"glm")
	figure2$MutType <- phenoname_all[i]
	return(data.frame(figure2))
},mc.cores=40))
write.csv(figure2_alltype,paste0(path,"/tableS3.csv"),row.names=F) #Figure 2B

## Figure 2C

model <- glm(DNM_Ms_0.35 ~ rcs(age_mat, quantile(age_mat, c(0.5,0.6,0.75))) + age_pat,data=dnm_final_foranalysis_singleton)
model_summary <- summary(model)

#tableS5
figure2c_data <- dnm_final_foranalysis_singleton

model1 <- lm(DNM_Ms_0.35~age_mat+age_pat,figure2c_data)
model2 <- lm(DNM_Ms_0.35~poly(age_mat,2)+age_pat,figure2c_data)
model3 <- lm(DNM_Ms_0.35~poly(age_mat,3)+age_pat,figure2c_data)

model4 <- lm(DNM_Ps_0.35~age_pat+age_mat,figure2c_data)
model5 <- lm(DNM_Ps_0.35~poly(age_pat,2)+age_mat,figure2c_data)
model6 <- lm(DNM_Ps_0.35~poly(age_pat,3)+age_mat,figure2c_data)

tableS5 <- data.frame(Group=c(rep("Maternal age",3),rep("Paternal age",3)),Modeltype=rep(c("Linear regression","Quadratic polynomial model","Cubic polynomial model"),2),Adjusted_R2=c(summary(model1)$adj.r.squared,summary(model2)$adj.r.squared,summary(model3)$adj.r.squared,summary(model4)$adj.r.squared,summary(model5)$adj.r.squared,summary(model6)$adj.r.squared),AIC=c(AIC(model1),AIC(model2),AIC(model3),AIC(model4),AIC(model5),AIC(model6)),P1=c("-",anova(model1,model2)$Pr[2],anova(model1,model3)$Pr[2],"-",anova(model4,model5)$Pr[2],anova(model4,model6)$Pr[2]),P2=c("-","-",anova(model2,model3)$Pr[2],"-","-",anova(model5,model6)$Pr[2]))
write.csv(tableS5,paste0(path,"/tableS5.csv"),row.names=F)

##
model1 <- lm(DNM_Ms_0.35~age_mat+age_pat,subset(dnm_final_foranalysis_singleton,age_mat<28.92))
model2 <- lm(DNM_Ms_0.35~age_mat+age_pat,subset(dnm_final_foranalysis_singleton,age_mat>31.81))

model1 <- glm(DNM_Ms_0.35 ~ rcs(age_mat, quantile(age_mat, c(0.5,0.6,0.75))) + age_pat,data=subset(dnm_final_foranalysis_singleton,age_mat<28.92))
model_summary <- summary(model1)

model2 <- glm(DNM_Ms_0.35 ~ rcs(age_mat, quantile(age_mat, c(0.5,0.6,0.75))) + age_pat,data=subset(dnm_final_foranalysis_singleton,age_mat>31.81))
model_summary <- summary(model2)

## Figure 2D

MO <- subset(dnm_final_foranalysis_singleton,age_mat>31.81)
MY <- subset(dnm_final_foranalysis_singleton,age_mat<28.92)

DNM_mut_all_final$indi <- ifelse(DNM_mut_all_final$VCFID%in%MO$IID,"MO",ifelse(DNM_mut_all_final$VCFID%in%MY$IID,"MY",""))

DNM_MO <- subset(DNM_mut_all_final,indi=="MO")
DNM_MY <- subset(DNM_mut_all_final,indi=="MY")

fisher_result2 <- data.frame(type=mutation_type2,num_type_MO=NA,num_all_MO=NA,num_type_MY=NA,num_all_MY=NA,ER=NA,LC=NA,UC=NA,P=NA)

for (i in 1:nrow(fisher_result2)) {
	fisher_result2$num_type_MO[i] <- nrow(subset(DNM_MO,tricontext==mutation_type2[i]))
	fisher_result2$num_all_MO[i] <- nrow(DNM_MO)
	fisher_result2$num_type_MY[i] <- nrow(subset(DNM_MY,tricontext==mutation_type2[i]))
	fisher_result2$num_all_MY[i] <- nrow(DNM_MY)
	test <- fisher.test(matrix(c(fisher_result2$num_type_MO[i],fisher_result2$num_all_MO[i]-fisher_result2$num_type_MO[i],fisher_result2$num_type_MY[i],fisher_result2$num_all_MY[i]-fisher_result2$num_type_MY[i]),ncol=2))
	fisher_result2$ER[i] <- test$estimate
	fisher_result2$LC[i] <- test$conf.int[1]
	fisher_result2$UC[i] <- test$conf.int[2]
	fisher_result2$P[i] <- test$p.value
}
subset(fisher_result2,grepl("C>G",type))
write.csv(fisher_result2,paste0(path,"/Figure2D.csv"),row.names=F)


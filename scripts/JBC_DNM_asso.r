options(stringsAsFactors=F)

source("JBC_DNM_func.r")

library(data.table)
library(parallel)

# ====================================================================
# Notes:
# 1. Linear regression results are reported as Beta (95% CI)
# 2. Logistic regression results are reported as Odds Ratio (95% CI)
# 3. All models are adjusted for predefined covariates
# 4. Reference category = first level of factor variables
# ====================================================================

# ====================================================================
# Part 1. Association between parental age and parental-origin DNMs
# ====================================================================

# Covariates

covar_list = c("age_mat","age_pat")
covar_class = rep("c",2)
covar_name = c("Maternal age","Paternal age")

# Outcomes (paternal-origin DNMs)

pheno_Ps = c("dSNVP","dIndelP","DNMP",paste0("dSNVP_",c("TG","TC","TA","CA","CG","CT_CpG","CT_nCpG")))
pheno_Ms = c("dSNVM","dIndelM","DNMM",paste0("dSNVM_",c("TG","TC","TA","CA","CG","CT_CpG","CT_nCpG")))

rbindlist(lapply(pheno_Ps, function(p) {
    res_art <- generate_multi_linear(dnm_jbc_input,covar_list,covar_class,covar_name,p)
    res_art$pheno <- p
    res_art[2,] # Row 2 corresponds to "Paternal age"
}),fill = TRUE)

rbindlist(lapply(pheno_Ms, function(p) {
    res_art <- generate_multi_linear(dnm_jbc_input,covar_list,covar_class,covar_name,p)
    res_art$pheno <- p
    res_art[1,] # Row 1 corresponds to "Maternal  age"
}),fill = TRUE)

# =========================================================
# Maternal / paternal age acceleration (non-linearity test)
# =========================================================

# Compare linear vs polynomial terms

model1 <- lm(dSNVM~age_mat+age_pat,dnm_jbc_input)
model2 <- lm(dSNVM~poly(age_mat,2)+age_pat,dnm_jbc_input)
model3 <- lm(dSNVM~poly(age_mat,3)+age_pat,dnm_jbc_input)

model4 <- lm(dSNVP~age_pat+age_mat,dnm_jbc_input)
model5 <- lm(dSNVP~poly(age_pat,2)+age_mat,dnm_jbc_input)
model6 <- lm(dSNVP~poly(age_pat,3)+age_mat,dnm_jbc_input)

anova(model1,model2)
anova(model4,model5)

# ====================================================================
# Part 2. Association between ART and DNMs
# ====================================================================

covar_list = c("age_pat","age_mat","Conception","gender")
covar_class = c(rep("c",2),rep("f",2))
covar_name = c("Paternal age","Maternal age","Conception type","Gender")

phenoname = c("dSNV","dSNVP","dSNVM","dIndel","dIndelP","dIndelM","DNM","DNMP","DNMM")

rbindlist(lapply(phenoname, function(p) {
    res_art <- generate_multi_linear(dnm_jbc_input,covar_list,covar_class,covar_name,p)
    res_art$pheno <- p
    res_art[4,] # Row 4 corresponds to "Conception type"
}),fill = TRUE)

# DDR sensitivity analysis

covar_list = c("age_pat","age_mat","Conception","gender","DDR_father","DDR_mother")
covar_class = c(rep("c",2),rep("f",2),rep("c",2))
covar_name = c("Paternal age","Maternal age","Conception type","Gender","DDR_mut_father","DDR_mut_mother")

phenoname=c("dSNVP","dSNVM")

rbindlist(lapply(phenoname, function(p) {
    res_art <- generate_multi_linear(dnm_jbc_input,covar_list,covar_class,covar_name,p)
    res_art$pheno <- p
    res_art[4,]
}),fill = TRUE)

# ====================================================================
# Part 3. Association between ART precedures and parental DNMs
# ====================================================================

# Subset ART conceptions
dnm_jbc_input_ART <- subset(dnm_jbc_input,Conception=="ART")

# Univariable analysis

vars <- c("ART_type","ART_stage","Embryo","OI_regimen","Trigger")

covar_list  <- c("age_pat","age_mat","gender")
covar_class <- c("c","c","f")
covar_name  <- c("Paternal age","Maternal age","Gender")

ART_procedure_uni_Ps_res <- rbindlist(lapply(vars, function(v) {
	res <- generate_multi_linear(dnm_jbc_input_ART, c(covar_list, v), c(covar_class, "f"), c(covar_name, v), "dSNVP")
	res
}), fill = TRUE)  

ART_procedure_uni_Ms_res <- rbindlist(lapply(vars, function(v) {
	res <- generate_multi_linear(dnm_jbc_input_ART, c(covar_list, v), c(covar_class, "f"), c(covar_name, v), "dSNVM")
	res
}), fill = TRUE)  

# infertility diagnosis sensitivity analysis

dnm_jbc_input_ART$Female_diag <- ifelse(dnm_jbc_input_ART$Female_ovulation_failure=="Yes" | dnm_jbc_input_ART$Female_pelvic=="Yes" | dnm_jbc_input_ART$Female_uterine=="Yes","Yes","No")
dnm_jbc_input_ART$Male_diag <- ifelse(dnm_jbc_input_ART$Male_oligothenospermia=="Yes"  | dnm_jbc_input_ART$Male_azoospermia=="Yes","Yes","No")

covar_list  <- c("age_pat","age_mat","gender")
covar_class <- c("c","c","f")
covar_name  <- c("Paternal age","Maternal age","Gender")

var_Ps <- c("Male_diag","Male_oligothenospermia","Male_azoospermia")
var_Ms <- c("Female_diag","Female_ovulation_failure","Female_pelvic","Female_uterine")

ART_infertility_uni_Ps_res <- rbindlist(lapply(var_Ps, function(v) {
	res <- generate_multi_linear(dnm_jbc_input_ART, c(covar_list, v), c(covar_class, "f"), c(covar_name, v), "dSNVP")
	res
}), fill = TRUE)  

ART_infertility_uni_Ms_res <- rbindlist(lapply(var_Ms, function(v) {
	res <- generate_multi_linear(dnm_jbc_input_ART, c(covar_list, v), c(covar_class, "f"), c(covar_name, v), "dSNVM")
	res
}), fill = TRUE)  

# multivariable analysis

covar_list=c("age_mat","age_pat","gender","ART_type","ART_stage","Embryo","OI_regimen","Trigger")
covar_class=c("c","c",rep("f",6))
covar_name=c("Maternal age","Paternal age","Gender","Fertilization method","Embryo transfer stage","Embryo transfer cycle","Regimen of ovulation induction","Trigger strategies")

generate_multi_linear(dnm_jbc_input_ART,covar_list,covar_class,covar_name,"dSNVP")
generate_multi_linear(dnm_jbc_input_ART,covar_list,covar_class,covar_name,"dSNVM")

# stratification analysis for paternal DNMs

dnm_jbc_input$ART_type_GP <- ifelse(dnm_jbc_input$Conception=="SP","NC",ifelse(dnm_jbc_input$Conception=="ART" & dnm_jbc_input$ART_type=="ICSI","ICSI",ifelse(dnm_jbc_input$Conception=="ART" & dnm_jbc_input$ART_type=="IVF","IVF","Unknown")))
dnm_jbc_input$ART_type_GP <- factor(dnm_jbc_input$ART_type_GP,levels=c("NC","IVF","ICSI","Unknown"))

# All conceptions
covar_list=c("age_mat","age_pat","gender","ART_type_GP")
covar_class=c("c","c",rep("f",2))
covar_name=c("Maternal age","Paternal age","Gender","Fertilization method")

generate_multi_linear(dnm_jbc_input,covar_list,covar_class,covar_name,"dSNVP")

# ART conceptions
covar_list=c("age_mat","age_pat","gender","ART_type","ART_stage","Embryo","OI_regimen","Trigger")
covar_class=c("c","c",rep("f",6))
covar_name=c("Maternal age","Paternal age","Gender","Fertilization method","Embryo transfer stage","Embryo transfer cycle","Regimen of ovulation induction","Trigger strategies")

generate_multi_linear(subset(dnm_jbc_input_ART,Male_oligothenospermia=="Yes"),covar_list,covar_class,covar_name,"dSNVP")[5:6,]
generate_multi_linear(subset(dnm_jbc_input_ART,Male_oligothenospermia=="No"),covar_list,covar_class,covar_name,"dSNVP")[5:6,]

# stratification analysis for maternal DNMs

dnm_jbc_input$GnRHant_doseGP_all <- ifelse(dnm_jbc_input$Conception=="SP","NC",dnm_jbc_input$GnRHant_doseGP)
dnm_jbc_input$GnRHant_doseGP_all <- factor(dnm_jbc_input$GnRHant_doseGP_all,levels=c("NC","0","Half_dose","Full_dose","High_dose"))

# All conceptions
covar_list=c("age_mat","age_pat","gender","GnRHant_doseGP_all")
covar_class=c("c","c",rep("f",2))
covar_name=c("Maternal age","Paternal age","Gender","GnRHant_dose")

generate_multi_linear(dnm_jbc_input,covar_list,covar_class,covar_name,"dSNVM")

# ART conceptions
covar_list=c("age_mat","age_pat","gender","ART_type","ART_stage","Embryo","GnRHant_doseGP_all","Trigger")
covar_class=c("c","c",rep("f",6))
covar_name=c("Maternal age","Paternal age","Gender","Fertilization method","Embryo transfer stage","Embryo transfer cycle","GnRHant_doseGP_all","Trigger strategies")

generate_multi_linear(subset(dnm_jbc_input,Conception=="ART"),covar_list,covar_class,covar_name,"dSNVM")

# =========================================================
# Part 4. Association between Parental DNMs and offspring outcomes
# =========================================================

dnm_jbc_input$age_ave <- (dnm_jbc_input$age_mat+dnm_jbc_input$age_pat)/2
dnm_jbc_input$PTB_rg <- ifelse(dnm_jbc_input$PTB=="Yes",1,ifelse(dnm_jbc_input$PTB=="No",0,dnm_jbc_input$PTB))
dnm_jbc_input$LBW_rg <- ifelse(dnm_jbc_input$LBW=="Yes",1,ifelse(dnm_jbc_input$LBW=="No",0,dnm_jbc_input$LBW))

# average parental age

covar_list = c("age_ave","BMI_Mother","Conception","gender","Income_annual","Education","Residence")
covar_class=c(rep("c",2),rep("f",5))
covar_name = c("Parental age","BMI Mother","Conception","Gender","Annual income","Education","Residence")

phenoname_binary <- c("PTB_rg","LBW_rg")

rbindlist(lapply(phenoname_binary, function(p) {
	res <- generate_multi_logit(dnm_jbc_input, covar_list, covar_class, covar_name, p)[1,] # Row 1 corresponds to average parental age
	res$pheno <- p
	res
}), fill = TRUE)

phenoname_continous <- c("GW","bw_kg")

rbindlist(lapply(phenoname_continous, function(p) {
	res <- generate_multi_linear_digit3(dnm_jbc_input, covar_list, covar_class, covar_name, p)[1,] # Row 1 corresponds to average parental age
	res$pheno <- p
	res
}), fill = TRUE)

# ART

covar_list = c("age_pat","age_mat","BMI_Mother","Conception","gender","Income_annual","Education","Residence")
covar_class=c(rep("c",3),rep("f",5))
covar_name = c("Paternal age","Maternal age","BMI Mother","Conception","Gender","Annual income","Education","Residence")

phenoname_binary <- c("PTB_rg","LBW_rg")

rbindlist(lapply(phenoname_binary, function(p) {
	res <- generate_multi_logit(dnm_jbc_input, covar_list, covar_class, covar_name, p)[5,] # Row 5 corresponds to Conception type
	res$pheno <- p
	res
}), fill = TRUE)

phenoname_continous <- c("GW","bw_kg")

rbindlist(lapply(phenoname_continous, function(p) {
	res <- generate_multi_linear_digit3(dnm_jbc_input, covar_list, covar_class, covar_name, p)[5,] # Row 5 corresponds to Conception type
	res$pheno <- p
	res
}), fill = TRUE)

# parental DNMs

vars <- c("dSNVP","dSNVM")
covar_list = c("age_pat","age_mat","BMI_Mother","Conception","gender","Income_annual","Education","Residence")
covar_class=c(rep("c",3),rep("f",5))
covar_name = c("Paternal age","Maternal age","BMI Mother","Conception","Gender","Annual income","Education","Residence")

phenoname_binary <- c("PTB_rg","LBW_rg")

rbindlist(lapply(vars, function(v) {
	rbindlist(lapply(phenoname_binary, function(p) {
		res <- generate_multi_logit(dnm_jbc_input, c(covar_list, v), c(covar_class, "c"), c(covar_name, v), p)[18,] # Row 18 corresponds to parental DNMs
		res$pheno <- p
		res$exposure <- v
		res
	}), fill = TRUE)
}), fill = TRUE)

phenoname_continous <- c("GW","bw_kg")

rbindlist(lapply(vars, function(v) {
	rbindlist(lapply(phenoname_continous, function(p) {
		res <- generate_multi_linear_digit3(dnm_jbc_input, c(covar_list, v), c(covar_class, "c"), c(covar_name, v), p)[18,] # Row 18 corresponds to parental DNMs
		res$pheno <- p
		res$exposure <- v
		res
	}), fill = TRUE)
}), fill = TRUE)

## mediation analysis

set.seed(1234)
model.m <- lm(dSNVP~Conception+age_ave+gender+BMI_Mother+Income_annual+Education+Residence,dnm_jbc_input)
model.y <- lm(GW~dSNVP+Conception+age_ave+gender+BMI_Mother+Income_annual+Education+Residence,dnm_jbc_input)
mediate_GW_ave_Ps <- mediate(model.m,model.y,sims=1000,boot=T,treat="age_ave",mediator="dSNVP")

set.seed(1234)
model.m <- lm(dSNVP~Conception+age_mat+age_pat+gender+BMI_Mother+Income_annual+Education+Residence,subset(dnm_jbc_input,Conception=="SP" | ART_type=="ICSI"))
model.y <- lm(GW~dSNVP+Conception+age_mat+age_pat+gender+BMI_Mother+Income_annual+Education+Residence,subset(dnm_jbc_input,Conception=="SP" | ART_type=="ICSI"))
mediate_P_ICSI_GW <- mediate(model.m,model.y,sims=1000,boot=T,treat="Conception",mediator="dSNVP")

# =========================================================
# Part 5. Association between ART and EPZMs
# =========================================================

covar_list = c("age_pat","age_mat","Conception","gender")
covar_class = c(rep("c",2),rep("f",2))
covar_name = c("Paternal age","Maternal age","Conception type","Gender")

phenoname = c("dSNVmo","dIndelmo","DNMmo")

rbindlist(lapply(phenoname, function(p) {
    res_art <- generate_multi_linear(dnm_jbc_input,covar_list,covar_class,covar_name,p)
    res_art$pheno <- p
    res_art[4,]
}),fill = TRUE)

# univariable

vars <- c("ART_type","ART_stage","Embryo","OI_regimen","Trigger")

covar_list  <- c("age_pat","age_mat","gender")
covar_class <- c("c","c","f")
covar_name  <- c("Paternal age","Maternal age","Gender")

rbindlist(lapply(vars, function(v) {
	res <- generate_multi_linear(dnm_jbc_input_ART, c(covar_list, v), c(covar_class, "f"), c(covar_name, v), "dSNVmo")
	res
}), fill = TRUE)  

# multivariable

covar_list=c("age_mat","age_pat","gender","ART_type","ART_stage","Embryo","OI_regimen","Trigger")
covar_class=c("c","c",rep("f",6))
covar_name=c("Maternal age","Paternal age","Gender","Fertilization method","Embryo transfer stage","Embryo transfer cycle","Regimen of ovulation induction","Trigger strategies")

generate_multi_linear(dnm_jbc_input_ART,covar_list,covar_class,covar_name,"dSNVmo")

# stratification analysis for EPZMs

dnm_jbc_input$ART_stage_GP <- ifelse(dnm_jbc_input$Conception=="SP","NC",ifelse(dnm_jbc_input$Conception=="ART" & dnm_jbc_input$Embryo=="Fresh" & dnm_jbc_input$ART_stage=="Cleavage","Fresh-Cleavage",ifelse(dnm_jbc_input$Conception=="ART" & dnm_jbc_input$Embryo=="Fresh" & dnm_jbc_input$ART_stage=="Blastocyst","Fresh-Blastocyst",ifelse(dnm_jbc_input$Conception=="ART" & dnm_jbc_input$Embryo=="Frozen" & dnm_jbc_input$ART_stage=="Cleavage","Frozen-Cleavage",ifelse(dnm_jbc_input$Conception=="ART" & dnm_jbc_input$Embryo=="Frozen" & dnm_jbc_input$ART_stage=="Blastocyst","Frozen-Blastocyst",NA)))))
dnm_jbc_input$ART_stage_GP <- factor(dnm_jbc_input$ART_stage_GP,levels=c("NC","Fresh-Cleavage","Fresh-Blastocyst","Frozen-Cleavage","Frozen-Blastocyst"))

# All conceptions
covar_list=c("age_mat","age_pat","gender","ART_stage_GP")
covar_class=c("c","c",rep("f",2))
covar_name=c("Maternal age","Paternal age","Gender","Embryo")

generate_multi_linear(dnm_jbc_input,covar_list,covar_class,covar_name,"dSNVmo")

# ART conceptions
covar_list=c("age_mat","age_pat","gender","ART_type","ART_stage_GP","OI_regimen","Trigger")
covar_class=c("c","c",rep("f",5))
covar_name=c("Maternal age","Paternal age","Gender","Fertilization method","Embryo transfer stage","Regimen of ovulation induction","Trigger strategies")

generate_multi_linear(subset(dnm_jbc_input,Conception=="ART"),covar_list,covar_class,covar_name,"dSNVmo")





# =========================================================
# Figure 1C. Comparison of trinucleotide context
# =========================================================

# Input:
# DNM_info:
# Event-level de novo mutation dataset containing one record per DNM,
# including parental origin assignment (paternal, maternal, mosaic),
# mutation class, and trinucleotide context.
#
# Note:
# Event-level DNM annotations used for this analysis are derived from
# the controlled-access DNM VCF files (GVM001081).
#
# Analysis:
# Mutation proportions were compared among paternal DNMs, maternal DNMs, 
# and EPZMs using Fisher's exact tests.

# =========================================================
# Prepare datasets
# =========================================================

DNM_list <- list(
    FPOG      = subset(DNM_info, Source %in% c("paternal","maternal")),
    EPZM       = subset(DNM_info, Source == "mosaic"),
    Paternal  = subset(DNM_info, Source == "paternal"),
    Maternal  = subset(DNM_info, Source == "maternal")
)

tricontext_type <- sort(unique(DNM_mut_all_final$tricontext))

# =========================================================
# Fisher's exact tests
# =========================================================

figure1c_result <- data.frame(
    muttype = tricontext_type,
    ER1 = NA, P1 = NA,
    ER2 = NA, P2 = NA,
    prop_P = NA,
    prop_M = NA,
    prop_EPZM = NA
)

for(i in seq_along(tricontext_type)){

    mt <- tricontext_type[i]

    ## FPOG vs EPZM
    test1 <- fisher.test(matrix(c(
        sum(DNM_list$EPZM$tricontext==mt),
        sum(DNM_list$EPZM$tricontext!=mt),
        sum(DNM_list$FPOG$tricontext==mt),
        sum(DNM_list$FPOG$tricontext!=mt)
    ),2))

    ## paternal vs maternal
    test2 <- fisher.test(matrix(c(
        sum(DNM_list$Paternal$tricontext==mt),
        sum(DNM_list$Paternal$tricontext!=mt),
        sum(DNM_list$Maternal$tricontext==mt),
        sum(DNM_list$Maternal$tricontext!=mt)
    ),2))

    figure1c_result[i,c("ER1","P1")] <- c(test1$estimate,test1$p.value)
    figure1c_result[i,c("ER2","P2")] <- c(test2$estimate,test2$p.value)

    figure1c_result$prop_P[i]    <- mean(DNM_list$Paternal$tricontext==mt)
    figure1c_result$prop_M[i]    <- mean(DNM_list$Maternal$tricontext==mt)
    figure1c_result$prop_EPZM[i] <- mean(DNM_list$EPZM$tricontext==mt)
}

figure1c_result$padj1 <- p.adjust(figure1c_result$P1,"fdr")
figure1c_result$padj2 <- p.adjust(figure1c_result$P2,"fdr")

# =========================================================
# Prepare data for plotting
# =========================================================

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

# =========================================================
# Generate Figure 
# =========================================================
library(cowplot)

color = c("#FEE07F","#E42927","#ECC7C5","#A2CF63","#507CB7","#91C0DB")
figure1c_d2$Source <- factor(figure1c_d2$Source,levels=c("Paternal dnSNV","Maternal dnSNV","EPZM"))

g1_up <- ggplot(figure1c_d2, aes(basenear, prop,fill=std.mutcat))+
		geom_bar(stat = "identity", position=position_dodge(),width=0.7)+
		theme_bw()+theme_classic()+
		theme(axis.text.y=element_text(size=10,color="#030303",face="bold"),axis.title.y=element_text(size=12,face="bold"),axis.text.x=element_text(color="black",size=10,angle=90, hjust = 1),axis.title.x=element_text(face="bold",color="black",size=12),panel.grid.major=element_line(colour=NA))+ 
		theme(legend.position = "none")+
		theme(strip.text.x = element_text(face = "bold", color = "white", hjust = 0, size = 15), strip.background.x = element_rect(fill = color, linetype = "solid",color = "white"), strip.background.y = element_rect(color = "white"))+
		xlab(NULL)+ylab("Mutation Proportion")+ 
		scale_fill_manual(values = color,breaks=c("C>A","C>G","C>T","T>A","T>C","T>G"),labels=c("C>A","C>G","C>T","T>A","T>C","T>G"))+facet_grid(Source~std.mutcat)
g1_down <- ggplot(figure1c_d1, aes(x=basenear, y=Source, fill=log2(ER+1))) + geom_tile(color="white") + 
		theme_bw() + 
		theme(axis.text.y=element_text(color="black",size=10,hjust = 1),axis.title.y=element_text(size=12,face="bold"),axis.text.x=element_text(color="black",size=10,angle=90, hjust = 1),axis.title.x=element_text(face="bold",color="black",size=12),panel.grid.major=element_line(colour=NA)) +
		xlab("") + ylab("Mutation Proportion") +
		geom_text(label=figure1c_d1$anno, size=2) + 
		scale_fill_distiller(palette = "Spectral")+facet_grid(.~std.mutcat)

plot_grid(g1_up, g1_down,ncol=1, align='v', axis="lr",rel_heights=c(2,1))



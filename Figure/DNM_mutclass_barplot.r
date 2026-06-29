# =========================================================
# Figure 1B. Comparison of mutational spectra
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
# Figure 1D was generated using the same strategy as Figure 1B.

# =========================================================
# Prepare datasets
# =========================================================

DNM_list <- list(
    FPOG      = subset(DNM_info, Source %in% c("paternal","maternal")),
    EPZM       = subset(DNM_info, Source == "mosaic"),
    Paternal  = subset(DNM_info, Source == "paternal"),
    Maternal  = subset(DNM_info, Source == "maternal")
)

mutation_type <- c("C>A","C>G","C>T(CpG)","C>T(Other)","T>G","T>C","T>A")

# =========================================================
# Fisher's exact tests
# =========================================================

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

# =========================================================
# Prepare data for plotting
# =========================================================
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

# =========================================================
# Generate Figure 
# =========================================================

figure1b_data$muttype <- factor(figure1b_data$muttype,levels=c("C>A","C>G","C>T(CpG)","C>T(Other)","T>A","T>C","T>G"))
figure1b_data$Source <- factor(figure1b_data$Source,levels=c("Maternal dnSNV","Paternal dnSNV","EPZM"))

f <- ggplot(figure1b_data, aes(muttype, prop, fill=Source)) +
    geom_bar(stat = "identity", position=position_dodge(),color="grey",width=0.7) + 
    theme_bw() +
    theme_classic() +
    theme(
        axis.text.y = element_text(size=12,color="#030303",face="bold"),
        axis.title.y = element_text(size=12,face="bold"),
        axis.text.x=element_text(face="bold",color="black",size=12),
        axis.title.x=element_text(face="bold",color="black",size=12),
        panel.grid.major=element_line(colour=NA),
        legend.title = element_blank()
	) +
    xlab(NULL) +
    ylab("Mutation Proportion") + 
    scale_fill_manual(
        values=c("#CC1B21","#0066BD","#FBB255"),
        breaks=c("Maternal dnSNV","Paternal dnSNV","EPZM"),
        labels=c("Maternal dnSNV","Paternal dnSNV","EPZM")
    )

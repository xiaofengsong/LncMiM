#!/usr/bin/env Rscript

library(parallel)
source("miMod.R")
source("lncMod.R")
source("pcMod.R")

args <- commandArgs(T)

miRNA_Exp <- read.table(args[1])
rownames(miRNA_Exp) <- miRNA_Exp[,1]
miRNA_Exp <- miRNA_Exp[,-1]
lncRNA_Exp <- read.table(args[2])
rownames(lncRNA_Exp) <- lncRNA_Exp[,1]
lncRNA_Exp <- lncRNA_Exp[,-1]
pc_Exp <- read.table(args[3])
rownames(pc_Exp) <- pc_Exp[,1]
pc_Exp <- pc_Exp[,-1]

tpt <- read.table(args[4]) 


out <- args[5]
ws <- 0.25
mode <- 'miRNA'
method <- 'spearman';
if (length(args) > 5){
	mode <- args[6]
	if (length(args) > 6){
		ws <- args[7]
		if (length(args) > 7){
			method <- args[8]
		}
	}
}

sink(out)
triplets <- as.data.frame(matrix(nrow = 0, ncol = 11))
colnames(triplets) <- c("modulator", "effector", "target","low_correlation", "high_correlation","low_index","high_index", "fisrt","last","p-value", "fdr")
if (mode == 'lncRNA'){
	ids <- unique(tpt[,2])
	for (md in ids){
		cp <- tpt[tpt[,2] == md,c(1,3)]
		result <- lncMod(md = md,CP = cp, M.exp = lncRNA_Exp,E.exp=miRNA_Exp,T.exp = pc_Exp,ws,method=method)
		triplets <- rbind(triplets,result$triplets)
	}
} else if (mode == 'mRNA'){
	ids <- unique(tpt[,3])
        for (md in ids){
                cp <- tpt[tpt[,3] == md,c(1,2)]
		result <- pcMod(md = md,CP = cp, M.exp = pc_Exp,E.exp=miRNA_Exp,T.exp = lncRNA_Exp,ws,method=method)
		triplets <- rbind(triplets,result$triplets)
	 }

} else {
	ids <- unique(tpt[,1])
        for (md in ids){
                cp <- tpt[tpt[,1] == md,c(2,3)]
		result <- miMod(md = md,CP = cp, M.exp = miRNA_Exp,E.exp=lncRNA_Exp,T.exp = pc_Exp,ws,method=method)
		triplets <- rbind(triplets,result$triplets)	
        }
}
print(triplets)
sink()


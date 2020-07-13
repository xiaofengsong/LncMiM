miMod <- function (md, CP, M.exp, E.exp, T.exp, WS = 0.25, method = "spearman",iqr.filter = c(log2(1.5), log2(1.5), log2(1.5)), cor.MvsCP = c(-0.3,-0.3), cor.EvsT.dif = 0.3, cor.EvsT.abs = 0.3, CP.fc.filter = log2(1.5),CP.p.filter = 0.01, rand = 100, correction = "BH", cores = 1)
{
    iqr <- apply(M.exp, 1, IQR, na.rm = TRUE)
    data.M <- M.exp[iqr > iqr.filter[1], ]
    iqr <- apply(E.exp, 1, IQR, na.rm = TRUE)
    data.E <- E.exp[iqr > iqr.filter[2], ]
    iqr <- apply(T.exp, 1, IQR, na.rm = TRUE)
    data.T <- T.exp[iqr > iqr.filter[3], ]
    index.CP.E <- as.character(CP[, 1]) %in% rownames(data.E)
    index.CP.T <- as.character(CP[, 2]) %in% rownames(data.T)
    index.CP <- index.CP.E & index.CP.T
    CP <- CP[index.CP, ]

    if (dim(CP)[1] == 0) {
        tmptmpMCP <- list()
        tmptriplets <- as.data.frame(matrix(nrow = 0, ncol = 11))
        colnames(tmptriplets) <- c("modulator", "effector", "target","low_correlation", "high_correlation","low_index","high_index", "fisrt","last","p-value", "fdr")
        tmptmpMCP[["triplets"]] <- tmptriplets
        tmpfactor <- "aa"
        tmpfactorindex <- F
        tmptmpMCP[["initialnot"]] <- tmpfactor[tmpfactorindex]
        tmptmpMCP[["filterdnot"]] <- tmpfactor[tmpfactorindex]
        return(tmptmpMCP)
    }
    else {
        infun <- function(m) {
            if (m %in% rownames(data.M)) {
                MCP <- as.data.frame(matrix(nrow = 0, ncol = 11))
                colnames(MCP) <- c("modulator", "effector", "target","low_correlation", "high_correlation","low_index","high_index", "fisrt","last","p-value", "fdr")
                file <- paste("correlation_matrix_", m, ".Rdata",sep = "")
                if (file.exists(file)) {
                  load(file)
                }
                else {
		  datasorted <- sort(unlist(data.M[m,]))
		  AllCorr <- c()
		  WS <- as.numeric(WS)
		  WS_len <- round(length(datasorted)*WS)
		  for (i in 1:(length(datasorted)-WS_len+1)){
			  index <-  names(datasorted[i:(i+WS_len-1)]) 
	                  CORR <- apply(CP, 1, function(x) cor(t(data.T[as.character(x[2]),index]), t(data.E[as.character(x[1]),index]), use = "pairwise.complete.obs", method = method))
		          AllCorr <- rbind(AllCorr,CORR)
        	  }
		  LOW <- apply(AllCorr,2,min)
		  LOW_index <- apply(AllCorr,2,which.min)
	          HIGH <- apply(AllCorr,2,max)
		  HIGH_index <- apply(AllCorr,2,which.max)			
		  FIRST <- AllCorr[1,]
		  LAST <- AllCorr[dim(AllCorr)[1],]		

	        }
                delta.cor <- HIGH - LOW
                index.diff <- abs(delta.cor) > cor.EvsT.dif
                index.abs <- LOW >= cor.EvsT.abs | HIGH >= cor.EvsT.abs
                index <- index.diff & index.abs
                if (sum(index) != 0) {
                  tmp <- cbind(rep(m, sum(index)), CP[index,], LOW[index], HIGH[index],LOW_index[index],HIGH_index[index],FIRST[index],LAST[index])
                  ml <- dim(data.M)[2]
                  LOW.indexs <- WS_len
                  HIGH.indexs <- WS_len
                  deltar <- tmp[5] - tmp[4]
                  for (rr in 1:rand) {
                    LOW.rand.index <- sample(1:ml, LOW.indexs)
                    HIGH.rand.index <- sample(setdiff(1:ml, LOW.rand.index),HIGH.indexs)
                    LOW.rand <- apply(tmp[, 2:3], 1, function(x) cor(t(data.T[as.character(x[2]), LOW.rand.index]), t(data.E[as.character(x[1]),LOW.rand.index]), use = "pairwise.complete.obs", method = method))
                    HIGH.rand <- apply(tmp[, 2:3], 1, function(x) cor(t(data.T[as.character(x[2]),HIGH.rand.index]), t(data.E[as.character(x[1]),HIGH.rand.index]), use = "pairwise.complete.obs",method = method))
                    deltar <- cbind(deltar, HIGH.rand - LOW.rand)
                  }
                  p <- numeric()
                  for (r in 1:dim(deltar)[1]) {
                    if (deltar[r, 1] <= 0) {
                      tmpp <- sum(deltar[r, 2:(rand + 1)] < deltar[r,1])/rand
                    }
                    if (deltar[r, 1] > 0) {
                      tmpp <- sum(deltar[r, 2:(rand + 1)] > deltar[r,1])/rand
                    }
                    p <- c(p, tmpp)
                  }
                  fdr <- p.adjust(p, method = correction)
                  tmp <- cbind(tmp, p, fdr)
                  colnames(tmp) <- c("modulator", "effector", "target","low_correlation", "high_correlation","low_index","high_index", "fisrt","last","p-value", "fdr")
                  MCP <- rbind(MCP, tmp)
                }
            }
            else {
                if (m %in% rownames(M.exp)) {
                  MCP <- "filterdnot"
                }
                else {
                  MCP <- "initialnot"
                }
            }
            return(MCP)
        }
        tmpMCP <- mclapply(md, infun, mc.cores = cores)
        names(tmpMCP) <- md
        findex <- grep("filterdnot", tmpMCP)
        iindex <- grep("initialnot", tmpMCP)
        index <- setdiff(1:length(tmpMCP), c(findex, iindex))
        tmptmpMCP <- list()
        tmptri <- as.data.frame(matrix(nrow = 0, ncol = 11))
        colnames(tmptri) <- c("modulator", "effector", "target","low_correlation", "high_correlation","low_index","high_index", "fisrt","last","p-value", "fdr")
        for (tmpi in index) {
            tmptri <- rbind(tmptri, tmpMCP[[tmpi]])
        }
        tmptmpMCP[["triplets"]] <- tmptri
        tmptmpMCP[["initialnot"]] <- md[iindex]
        tmptmpMCP[["filterdnot"]] <- md[findex]
        return(tmptmpMCP)
    }
}


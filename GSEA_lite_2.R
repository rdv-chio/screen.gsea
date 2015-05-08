GSEA_lite_2 <- function(input.ds, gs.db, nperm = 1000, fdr.q.val.threshold = 0.1, output.directory = '', doc.string='GSEA.analysis', weighted.score.type=0) {
 
# Read input gene list (must be stored as data.frame)

	print("Reading input gene list")

	dataset <- read.table(input.ds,header=T,sep='\t')
	dataset <- dataset[order(dataset[[2]],decreasing=T),]
	dataset.genes <- as.character(dataset[,1])
	dataset.effects <- dataset[,2] # effect size must not equal zero
	dataset.genes.N <- length(dataset.genes)
	dataset.indexes <- c(1:dataset.genes.N)

# Read input gene set database (must be in .gmt format)

	print("Reading input gene set database")

	genesets <- readLines(gs.db)
	geneset.N <- length(genesets)
	geneset.sizes <- vector(length = geneset.N, mode = "numeric") 
	for (i in 1:geneset.N) {
		geneset.sizes[i] <- length(which(unlist(strsplit(genesets[[i]], "\t")) != 'null')) - 2
	}
	max.geneset.size <- max(geneset.sizes)      
	geneset.matrix <- matrix(rep("null", geneset.N*max.geneset.size), nrow=geneset.N, ncol= max.geneset.size)
	geneset.names <- vector(length = geneset.N, mode = "character")
	geneset.descriptions <- vector(length = geneset.N, mode = "character")
	for (i in 1:geneset.N) {
		geneset.line <- noquote(unlist(strsplit(genesets[[i]], "\t")))
		geneset.names[i] <- geneset.line[1] 
		geneset.descriptions[i] <- geneset.line[2] 
		geneset.genes <- vector(length = geneset.sizes[i], mode = "character")
		for (j in 1:geneset.sizes[i]) {
			geneset.genes[j] <- geneset.line[j + 2]
		}
		geneset.matrix[i,] <- c(geneset.genes, rep("null", max.geneset.size - geneset.sizes[i]))
	}

# To eliminate rows from geneset.matrix with zero components in gene list

	geneset.matrix.match=matrix(rep(NA,times=nrow(geneset.matrix)*ncol(geneset.matrix)),nrow(geneset.matrix),ncol(geneset.matrix))
	for (i in 1:nrow(geneset.matrix)) {
		for (j in 1:length(geneset.matrix[i,][geneset.matrix[i,]!='null'])) {
			geneset.matrix.match[i,j]=match(geneset.matrix[i,j],dataset.genes,0)
		}
	}
	match.sum=vector('numeric',length=nrow(geneset.matrix.match))
	for (i in 1:length(match.sum)) {
		match.sum[i]=sum(geneset.matrix.match[i,],na.rm=T)
	}
	geneset.good.matrix=matrix(rep(NA,times=nrow(geneset.matrix)*ncol(geneset.matrix)),nrow(geneset.matrix),ncol(geneset.matrix))
	row.names(geneset.good.matrix)=c(1:nrow(geneset.matrix))
	geneset.bad.matrix=matrix(rep(NA,times=nrow(geneset.matrix)*ncol(geneset.matrix)),nrow(geneset.matrix),ncol(geneset.matrix))
	for (i in 1:length(match.sum)) {
		if (match.sum[i]>0) {
			geneset.good.matrix[i,]=geneset.matrix[i,]
		} else {
			geneset.bad.matrix[i,]=geneset.matrix[i,]
		}
	}
	geneset.matrix=geneset.good.matrix[rowSums(is.na(geneset.good.matrix))!=ncol(geneset.good.matrix),]
	geneset.matrix.match.2=matrix(rep(NA,times=nrow(geneset.matrix)*ncol(geneset.matrix)),nrow(geneset.matrix),ncol(geneset.matrix))
	for (i in 1:nrow(geneset.matrix)) {
		for (j in 1:length(geneset.matrix[i,][geneset.matrix[i,]!='null'])) {
			geneset.matrix.match.2[i,j]=match(geneset.matrix[i,j],dataset.genes,0)
		}
	}
	geneset.hit.matrix=matrix(rep("null",times=nrow(geneset.matrix)*ncol(geneset.matrix)),nrow(geneset.matrix),ncol(geneset.matrix))
	geneset.hit.vector=vector('character',length=nrow(geneset.matrix))
	geneset.miss.matrix=matrix(rep("null",times=nrow(geneset.matrix)*ncol(geneset.matrix)),nrow(geneset.matrix),ncol(geneset.matrix))
	geneset.miss.vector=vector('character',length=nrow(geneset.matrix))
	geneset.hit.matrix[c(which(geneset.matrix.match.2>0))] = geneset.matrix[c(which(geneset.matrix.match.2>0))]
	geneset.miss.matrix[c(which(geneset.matrix.match.2==0))] = geneset.matrix[c(which(geneset.matrix.match.2==0))]
	for (i in 1:nrow(geneset.hit.matrix)) {
		geneset.hit.vector[i] = paste(geneset.hit.matrix[i,geneset.hit.matrix[i,]!='null'],collapse=";")
		geneset.miss.vector[i] = paste(geneset.miss.matrix[i,geneset.miss.matrix[i,]!='null'],collapse=";")
	}
	geneset.N=nrow(geneset.matrix)
	geneset.names.2=vector(length = geneset.N, mode = "character")
	geneset.descriptions.2=vector(length = geneset.N, mode = "character")
	geneset.sizes.2=vector(length = geneset.N, mode = "numeric")
	for (i in 1:geneset.N) {
		geneset.names.2[i]=geneset.names[as.numeric(row.names(geneset.matrix)[i])]
		geneset.descriptions.2[i]=geneset.descriptions[as.numeric(row.names(geneset.matrix)[i])]
		geneset.sizes.2[i]=geneset.sizes[as.numeric(row.names(geneset.matrix)[i])]
	}

#Compute observed enrichment scores for each gene set

	print("Computing enrichment scores for observed gene list")
	
	obs.ES <- vector(length = geneset.N, mode = "numeric")
	obs.arg.ES <- vector(length = geneset.N, mode = "numeric")
	obs.ES.norm <- vector(length = geneset.N, mode = "numeric")
	obs.RES <- matrix(nrow= geneset.N, ncol=dataset.genes.N)
	obs.indicator <- matrix(nrow= geneset.N, ncol=dataset.genes.N)
	
	for (i in 1:geneset.N) {
		geneset <- geneset.matrix[i,geneset.matrix[i,] != "null"]
		geneset.indexes <- vector(length=length(geneset), mode = "numeric")
		geneset.indexes <- match(geneset, dataset.genes)
		GSEA.results <- GSEA_ES_lite_2(gene.list=dataset.indexes, gene.set=geneset.indexes, weighted.score.type, correl.vector=dataset.effects)
		obs.ES[i] <- GSEA.results$ES
		obs.arg.ES[i] <- GSEA.results$arg.ES
		obs.RES[i,] <- GSEA.results$RES
		obs.indicator[i,] <- GSEA.results$indicator
	}

# Compute enrichment scores for random permutations

	print("Computing enrichment scores for randomly permuted gene lists")

	phi <- matrix(nrow = geneset.N, ncol = nperm)
	phi.norm <- matrix(nrow = geneset.N, ncol = nperm)
	obs.phi <- matrix(nrow = geneset.N, ncol = nperm)
	obs.phi.norm <- matrix(nrow = geneset.N, ncol = nperm)

	for (i in 1:geneset.N) {
		geneset <- geneset.matrix[i,geneset.matrix[i,] != "null"]
		geneset.indexes <- vector(length=length(geneset), mode = "numeric")	
		geneset.indexes <- match(geneset, dataset.genes)
		for (r in 1:nperm) {
			reshuffled.dataset.indexes <- sample(1:dataset.genes.N)
			GSEA.results <- GSEA_ES_lite_2(gene.list=reshuffled.dataset.indexes, gene.set=geneset.indexes, weighted.score.type, correl.vector=dataset.effects)   
			phi[i, r] <- GSEA.results$ES
		}
		GSEA.results <- GSEA_ES_lite_2(gene.list=dataset.indexes, gene.set=geneset.indexes, weighted.score.type, correl.vector=dataset.effects)
		obs.phi[i, 1] <- GSEA.results$ES
		for (r in 2:nperm) {
			obs.phi[i, r] <- obs.phi[i, 1]
		}
	}

# Normalize observed and permuted enrichment scores

	print("Normalizing enrichment scores")

	for (i in 1:geneset.N) {
		pos.phi <- NULL
		neg.phi <- NULL
		for (j in 1:nperm) {
			if (phi[i, j] >= 0) {
				pos.phi <- c(pos.phi, phi[i, j])
			} else {
				neg.phi <- c(neg.phi, phi[i, j])
			}
		}
		pos.m <- mean(as.numeric(pos.phi))
		neg.m <- mean(abs(as.numeric(neg.phi)))
		for (j in 1:nperm) {
			if (phi[i, j] >= 0) {
				phi.norm[i, j] <- phi[i, j]/pos.m 
			} else {
				phi.norm[i, j] <- phi[i, j]/neg.m
			}
		}
		for (j in 1:nperm) {
			if (obs.phi[i, j] >= 0) {
				obs.phi.norm[i, j] <- obs.phi[i, j]/pos.m
			} else {
				obs.phi.norm[i, j] <- obs.phi[i, j]/neg.m
			}
		}
		if (obs.ES[i] >= 0) {
			obs.ES.norm[i] <- obs.ES[i]/pos.m
		} else {
			obs.ES.norm[i] <- obs.ES[i]/neg.m
		}
	}

# Compute FDRs 

	print("Computing FDR q-values")

	phi.norm.mean <- vector(length=geneset.N, mode="numeric")
	obs.phi.norm.mean  <- vector(length=geneset.N, mode="numeric")
	obs.phi.mean  <- vector(length=geneset.N, mode="numeric")
	FDR.mean <- vector(length=geneset.N, mode="numeric")
	obs.ES.index <- order(obs.ES.norm, decreasing=T) # position of each normalized ES score in rank of ES scores, in decreasing order
	orig.index <- seq(1, geneset.N) # original geneset indexes
	orig.index <- orig.index[obs.ES.index] # original geneset indexes ranked in order of decreasing ES scores
	orig.index <- order(orig.index, decreasing=F) # rank of each geneset (in original order) according to decreasing ES score
	obs.ES.norm.sorted <- obs.ES.norm[obs.ES.index] # normalized ES values sorted in decreasing order
	geneset.names.sorted <- geneset.names.2[obs.ES.index] # geneset names sorted according to ES score rank, in decreasing order 

	for (k in 1:geneset.N) {
		ES.value <- obs.ES.norm.sorted[k]
		count.col <- vector(length=nperm, mode="numeric")
		obs.count.col <- vector(length=nperm, mode="numeric")
		for (i in 1:nperm) {
			phi.vec <- phi.norm[,i]
			obs.phi.vec <- obs.phi.norm[,i]
			if (ES.value >= 0) {
				count.col.norm <- sum(phi.vec >= 0) # add up positive ES scores for the permuted gene list across all genesets
				obs.count.col.norm <- sum(obs.phi.vec >= 0) # add up positive ES scores for the input gene list across all genesets
				count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec >= ES.value)/count.col.norm, 0) # report fraction of the sum of positive ES scores for the permuted gene list across all genesets that is larger than the ES score for the input gene list for geneset[k]
				obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 0) # report fraction of the sum of positive ES scores for the input gene list across all genesets that is larger than the ES score for the input gene list for geneset[k]
			} else {
				count.col.norm <- sum(phi.vec < 0)
				obs.count.col.norm <- sum(obs.phi.vec < 0)
				count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec <= ES.value)/count.col.norm, 0)
				obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 0)
			}
		}
		phi.norm.mean[k] <- mean(count.col) # for geneset[k] across all permutations, average the fraction of the sum of positive ES scores for the permuted gene list across all genesets that is larger than the ES score for the input gene list for geneset[k]
		obs.phi.norm.mean[k] <- mean(obs.count.col) # for geneset[k] across all permutations, average the fraction of the sum of positive ES scores for the input gene list across all genesets that is larger than the ES score for the input gene list for geneset[k]
		FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, phi.norm.mean[k]/obs.phi.norm.mean[k], 1) # for geneset[k], report ratio of phi.norm.mean[k] to obs.phi.norm.mean[k]
	}
	FDR.mean.sorted <- FDR.mean[orig.index]

# Produce report and running enrichment plot for each gene set passing the FDR q-value cut-off

	print("Producing result tables and plots")
	
	for (i in 1:geneset.N) {
		if (FDR.mean.sorted[i] <= fdr.q.val.threshold) {

		# Produce report

			kk <- 1
			gene.number <- vector(length = geneset.sizes.2[i], mode = "character")
			gene.name <- vector(length = geneset.sizes.2[i], mode = "character")
			gene.list.loc <- vector(length = geneset.sizes.2[i], mode = "numeric")
			core.enrichment <- vector(length = geneset.sizes.2[i], mode = "character")
			gene.effect <- vector(length = geneset.sizes.2[i], mode = "numeric")
			gene.RES <- vector(length = geneset.sizes.2[i], mode = "numeric")
			rank.list <- seq(1, dataset.genes.N)
			if (obs.ES[i] >= 0) {
				set.k <- seq(1, dataset.genes.N, 1)
				loc <- match(i, obs.ES.index)
			} else {
				set.k <- seq(dataset.genes.N, 1, -1)
				loc <- geneset.N - match(i, obs.ES.index) + 1
			}
			for (k in set.k) {
				if (obs.indicator[i, k] == 1) {
					gene.number[kk] <- kk
					gene.name[kk] <- dataset.genes[k]
					gene.list.loc[kk] <- k
					gene.effect[kk] <- dataset.effects[k]
					gene.RES[kk] <- signif(obs.RES[i, k], digits = 3)
					if (obs.ES[i] >= 0) {
						core.enrichment[kk] <- ifelse(gene.list.loc[kk] <= obs.arg.ES[i], "YES", "NO")
					} else {
						core.enrichment[kk] <- ifelse(gene.list.loc[kk] > obs.arg.ES[i], "YES", "NO")
					}
					kk <- kk + 1
				}
			}
			gene.report <- data.frame(cbind(gene.number, gene.name, gene.list.loc, gene.effect, gene.RES, core.enrichment))
			names(gene.report) <- c("#", "GENE", "LIST LOCATION", "EFFECT", "RES", "CORE_ENRICHMENT")
			gene.report <- gene.report[gene.report[[2]] != '',]
			if (output.directory != "")  {
				filename <- paste(output.directory, doc.string, ".", substring(geneset.names.2[i],1,30), ".report", ".", loc, ".txt", sep="", collapse="")
				write.table(gene.report, file = filename, quote=F, row.names=F, sep = "\t")
				if (.Platform$OS.type == "unix") {
					gs.filename <- paste(output.directory, doc.string, ".",substring(geneset.names.2[i],1,30), ".plot", ".", loc, ".pdf", sep="", collapse="")
					pdf(file=gs.filename, height = 6, width = 14)
				} else if (.Platform$OS.type == "windows") {
					gs.filename <- paste(output.directory, doc.string, ".", substring(geneset.names.2[i],1,30), ".plot", ".", loc, ".pdf", sep="", collapse="")
					pdf(file=gs.filename, height = 6, width = 14)
				}
			}

		# Produce running enrichment plot

			nf <- layout(matrix(c(1,2,3), 1, 3, byrow=T), 1, c(1, 1, 1), TRUE)
			ind <- 1:dataset.genes.N
			min.RES <- min(obs.RES[i,])
			max.RES <- max(obs.RES[i,])
			if (max.RES < 0.3) max.RES <- 0.3
			if (min.RES > -0.3) min.RES <- -0.3
			delta <- (max.RES - min.RES)*0.50
			min.plot <- min.RES - 2*delta
			max.plot <- max.RES
			max.corr <- max(dataset.effects)
			min.corr <- min(dataset.effects)
			obs.correl.vector.norm <- (dataset.effects - min.corr)/(max.corr - min.corr)*1.25*delta + min.plot
			zero.corr.line <- (- min.corr/(max.corr - min.corr))*1.25*delta + min.plot
			col <- ifelse(obs.ES[i] > 0, 2, 4)
			sub.string <- paste("Number of genes: ", dataset.genes.N, " (in list), ", geneset.sizes.2[i], " (in gene set)", sep = "", collapse="")
			main.string <- paste("Gene Set ", i, ":", geneset.names.2[i])
			plot(ind, obs.RES[i,], main = main.string, sub = sub.string, xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)", xlim=c(1, dataset.genes.N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, cex.main = 1, col = col)
			for (j in seq(1, dataset.genes.N, 20)) {
				lines(c(j, j), c(zero.corr.line, obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
			}
			lines(c(1, dataset.genes.N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
			lines(c(obs.arg.ES[i], obs.arg.ES[i]), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
			for (j in 1:dataset.genes.N) {
				if (obs.indicator[i, j] == 1) {
					lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
				}
			}
			lines(ind, obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
			lines(c(1, dataset.genes.N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
			leg.txt <- paste("Peak at ", obs.arg.ES[i], sep="", collapse="")
			text(x=obs.arg.ES[i], y=min.plot + 1.8*delta, labels=leg.txt, cex = 1.0)
			dev.off()
		}
	}

# Produce global results report

	obs.ES <- signif(obs.ES, digits=5)
	obs.ES.norm <- signif(obs.ES.norm, digits=5)
	FDR.mean.sorted <- signif(FDR.mean.sorted, digits=5)
	hit.count <- rowSums(obs.indicator)
	hit.percent <- hit.count/geneset.sizes.2

	global.report <- data.frame(cbind(geneset.names.2, geneset.descriptions.2, geneset.sizes.2, hit.count, hit.percent, geneset.hit.vector, geneset.miss.vector, obs.ES, obs.ES.norm, FDR.mean.sorted))
	names(global.report) <- c("GENESET", "DESCRIPTION", "SIZE", "HIT#", "HIT%", "HITS", "MISSES", "ES", "NES", "FDR Q-VAL")
	global.report <- global.report[order(global.report$"FDR Q-VAL"),]
	if (output.directory != "")  {
		global.filename <- paste(output.directory, doc.string, ".global.report.txt", sep="", collapse="")
		write.table(global.report, file = global.filename, quote=F, row.names=F, sep = "\t")
	}
}
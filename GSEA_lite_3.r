GSEA_lite_3 <- function(input.ds, gs.db, output.directory, nperm = 1000, weight=0) {
 
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
			geneset.matrix[i,j] <- geneset.line[j + 2]
		}
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
		match.sum[i]=length(geneset.matrix.match[i,c(which(geneset.matrix.match[i,]>0))])
	}
	geneset.good.matrix=matrix(rep(NA,times=nrow(geneset.matrix)*ncol(geneset.matrix)),nrow(geneset.matrix),ncol(geneset.matrix))
	row.names(geneset.good.matrix)=c(1:nrow(geneset.matrix))
	geneset.bad.matrix=matrix(rep(NA,times=nrow(geneset.matrix)*ncol(geneset.matrix)),nrow(geneset.matrix),ncol(geneset.matrix))
	for (i in 1:length(match.sum)) {
		if (match.sum[i]>3) {
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

#Compute enrichment scores for each gene set in original and permuted gene lists

	print("Computing enrichment scores")
	
	obs.ES <- vector(length = geneset.N, mode = "numeric")
	obs.med.ES <- vector(length = geneset.N, mode = "numeric")
	obs.arg.ES <- vector(length = geneset.N, mode = "numeric")
	obs.arg.med.ES <- vector(length = geneset.N, mode = "numeric")
	obs.RES <- matrix(nrow= geneset.N, ncol=dataset.genes.N)
	obs.indicator <- matrix(nrow= geneset.N, ncol=dataset.genes.N)
	phi <- matrix(nrow = geneset.N, ncol = nperm)
	phi.med <- matrix(nrow = geneset.N, ncol = nperm)
	
	for (i in 1:geneset.N) {
		geneset <- geneset.matrix[i,geneset.matrix[i,] != "null"]
		geneset.indexes <- vector(length=length(geneset), mode = "numeric")
		geneset.indexes <- match(geneset, dataset.genes)
		GSEA.results <- GSEA_ES_lite_3(gene.list=dataset.indexes, gene.set=geneset.indexes, weight, correl.vector=dataset.effects)
		obs.ES[i] <- GSEA.results$ES
		obs.med.ES[i] <- GSEA.results$med.ES
		obs.arg.ES[i] <- GSEA.results$arg.ES
		obs.arg.med.ES[i] <- GSEA.results$arg.med.ES
		obs.RES[i,] <- GSEA.results$RES
		obs.indicator[i,] <- GSEA.results$indicator
		for (r in 1:nperm) {
			reshuffled.dataset.indexes <- sample(1:dataset.genes.N)
			GSEA.results <- GSEA_ES_lite_3(gene.list=reshuffled.dataset.indexes, gene.set=geneset.indexes, weight, correl.vector=dataset.effects)   
			phi[i, r] <- GSEA.results$ES
			phi.med[i, r] <- GSEA.results$med.ES
		}
	}

# Normalize observed and permuted enrichment scores

	print("Normalizing enrichment scores")
	
	phi.norm <- matrix(nrow = geneset.N, ncol = nperm)
	phi.med.norm <- matrix(nrow = geneset.N, ncol = nperm)
	obs.ES.norm <- vector(length = geneset.N, mode = "numeric")
	obs.med.ES.norm <- vector(length = geneset.N, mode = "numeric")
	
	for (i in 1:geneset.N) {
		mean.phi <- mean(as.numeric(phi[i,]))
		mean.phi.med <- mean(as.numeric(phi.med[i,]))
		phi.norm[i,] <- phi[i,]/mean.phi
		obs.ES.norm[i] <- obs.ES[i]/mean.phi
		phi.med.norm[i,] <- phi.med[i,]/mean.phi.med
		obs.med.ES.norm[i] <- obs.med.ES[i]/mean.phi.med
	}

# Compute FDRs 

	print("Computing FDR q-values")

	FDR.ES <- vector(length=geneset.N, mode="numeric")
	FDR.med.ES <- vector(length=geneset.N, mode="numeric")

	for (i in 1:geneset.N) {
		FDR.ES[i] <- sum(phi.norm[i,] >= obs.ES.norm[i])/sum(phi.norm[i,] >= 0) # compute the fraction of ES scores larger than that observed for the original gene list across all permutations
		FDR.med.ES[i] <- sum(phi.med.norm[i,] >= obs.med.ES.norm[i])/sum(phi.med.norm[i,] >= 0) # compute the fraction of ES scores larger than that observed for the original gene list across all permutations
	}

# Produce report and running enrichment plot for each gene set passing the FDR q-value cut-off

	print("Producing result tables and plots")

	for (i in 1:geneset.N) {
	
		#Produce gene report
	
		gene.number <- vector(length = geneset.sizes.2[i], mode = "character")
		gene.name <- vector(length = geneset.sizes.2[i], mode = "character")
		gene.list.loc <- vector(length = geneset.sizes.2[i], mode = "numeric")
		gene.effect <- vector(length = geneset.sizes.2[i], mode = "numeric")
		gene.RES <- vector(length = geneset.sizes.2[i], mode = "numeric")
		core.enrichment <- vector(length = geneset.sizes.2[i], mode = "character")
		
		kk <- 1
		set.k <- seq(1, dataset.genes.N, 1)
		
		for (k in set.k) {
			if (obs.indicator[i, k] == 1) {
				gene.number[kk] <- kk
				gene.name[kk] <- dataset.genes[k]
				gene.list.loc[kk] <- k
				gene.effect[kk] <- dataset.effects[k]
				gene.RES[kk] <- signif(obs.RES[i, k], digits = 3)
				core.enrichment[kk] <- ifelse(gene.list.loc[kk] <= obs.arg.ES[i], "YES", "NO")
				kk <- kk + 1
			}
		}
		
		gene.report <- data.frame(cbind(gene.number, gene.name, gene.list.loc, gene.effect, gene.RES, core.enrichment))
		names(gene.report) <- c("#", "GENE", "LIST LOCATION", "EFFECT", "RES", "CORE_ENRICHMENT")
		gene.report <- gene.report[gene.report[[2]] != '',]
		filename <- paste(output.directory, substring(geneset.names.2[i],1,40), ".report", ".txt", sep="", collapse="")
		write.table(gene.report, file = filename, quote=F, row.names=F, sep = "\t")
		
		#Produce running enrichment plot
		
		pdf.filename <- paste(output.directory,substring(geneset.names.2[i],1,40), ".plot", ".pdf", sep="", collapse="")
		pdf(file=pdf.filename, height = 6, width = 14)
		
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
		sub.string <- paste("Number of genes: ", dataset.genes.N, " (in list), ", geneset.sizes.2[i], " (in gene set)", sep = "", collapse="")
		main.string <- geneset.names.2[i]
		plot(ind, obs.RES[i,], main = main.string, sub = sub.string, xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)", xlim=c(1, dataset.genes.N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, cex.main = 1, col = 2)
		for (j in seq(1, dataset.genes.N, 10)) {
			lines(c(j, j), c(zero.corr.line, obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
		}
		lines(c(1, dataset.genes.N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
		lines(c(obs.arg.ES[i], obs.arg.ES[i]), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 2) # max enrichment vertical line
		lines(c(obs.arg.med.ES[i], obs.arg.med.ES[i]), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 4) # median enrichment vertical line
		for (j in 1:dataset.genes.N) {
			if (obs.indicator[i, j] == 1) {
				lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
			}
		}
		lines(ind, obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
		lines(c(1, dataset.genes.N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
		leg.txt <- paste("Max at ", obs.arg.ES[i], sep="", collapse="")
		text(x=obs.arg.ES[i], y=min.plot + 2*delta, labels=leg.txt, cex = 1.0)
		leg.txt <- paste("Median at ", obs.arg.med.ES[i], sep="", collapse="")
		text(x=obs.arg.med.ES[i], y=min.plot + 2.2*delta, labels=leg.txt, cex = 1.0)
		dev.off()
	}
	
# Produce global results report

	obs.ES <- signif(obs.ES, digits=5)
	obs.med.ES <- signif(obs.med.ES, digits=5)
	obs.ES.norm <- signif(obs.ES.norm, digits=5)
	obs.med.ES.norm <- signif(obs.med.ES.norm, digits=5)
	FDR.ES <- signif(FDR.ES, digits=5)
	FDR.med.ES <- signif(FDR.med.ES, digits=5)
	hit.count <- rowSums(obs.indicator)
	hit.percent <- hit.count/geneset.sizes.2

	global.report <- data.frame(cbind(geneset.names.2, geneset.descriptions.2, geneset.sizes.2, geneset.hit.vector, obs.ES, obs.ES.norm, FDR.ES, obs.med.ES, obs.med.ES.norm, FDR.med.ES))
	names(global.report) <- c("GENESET", "DESCRIPTION", "SIZE", "HITS", "ES", "NES", "FDR", "MED ES", "MED NES", "MED FDR")
	global.report <- global.report[order(global.report$"FDR"),]
	if (output.directory != "")  {
		global.filename <- paste(output.directory, "global.report.txt", sep="", collapse="")
		write.table(global.report, file = global.filename, quote=F, row.names=F, sep = "\t")
	}
}
GSEA = function(input, annotation_file, geneset_database, output.directory, pval_threshold = 1, fdr_threshold = .1, absolute_effects = "TRUE", effect_type = "mean", ES_type = "med", geneset_size_threshold = 3, nperm = 1000, weight = 1) {

## input: Input filename. Columns of input should be: shrna, gene, effect, pval.
## annotation_file: Current gene annotations for each hairpin. Columns should be: shrna, gene.
## geneset_database: Genesets to be analyzed.
## output.directory: Directory where results should be stored.
## pval_threshold: Numeric 0-1. Hairpins with p-values greater than this threshold will be removed from the dataset.
## fdr_threshold: Numeric 0-1. Genesets with FDR values greater than this threshold will be excluded from summary plots.
## absolute_effects: Logical. If "TRUE", all effects will be converted to positive values. 
## effect_type: "mean" or "max". Determines whether mean or maximum effect should be used for each gene.
## ES_type: "med" or "max". Determines whether the ES for a geneset is determined by the median or the maximum of the RES.
## geneset.size.threshold: Numeric. Genesets with fewer hits in the gene list than this threshold will be excluded from the analysis.
## nperm: Numeric. Determines the number of permutations to be performed during the FDR calculation.
## weight: Numeric. This value is used to exponentially scale the effect sizes of each hairpin during ES calculation.

## Read dataset

	print("Curating gene list")

	data = read.table(input, header = T, sep = '\t', colClasses = c('character', 'character', 'numeric', 'numeric'))
	
## Update shRNA annotations
	
	annotations = read.table(annotation_file, header = T, sep = '\t', colClasses = c('character', 'character'))
	match.vector = match(annotations[[1]], data[[1]], nomatch = 0)
	annotations = annotations[match.vector > 0,]
	match.vector = match.vector[match.vector > 0]
	data.effects = data[match.vector,3]
	data.pvals = data[match.vector,4]
	data = data.frame(cbind(annotations, data.effects, data.pvals), stringsAsFactors = F)
	data = data[which(!is.na(data[[3]])),]
	colnames(data) = c('shRNA', 'Gene', 'Effect', 'P Value')
	
## Filter based on pval_threshold
	
	data = data[data[[4]] <= pval_threshold,]
	
## Calculate mean effect for each gene

	genes = unique(data[[2]])
	mean.effects = vector('numeric', length(genes))
	for(i in 1:length(genes)) {
		effects = data[data[[2]] == genes[i],3]
		mean.effects[i] = mean(effects)
	}
	if(absolute_effects == "TRUE") {
		mean.effects = abs(mean.effects)
	} else {
	}
		
## Calculate max effect for each gene (genes with both positive and negative effects are removed)

	both.indicator = vector('numeric', length(genes))
	max.effects = vector('numeric', length(genes))
	min.effects = vector('numeric', length(genes))
	for(i in 1:length(genes)) {
		effects = data[data[[2]] == genes[i],3]
		max.effects[i] = max(effects)
		min.effects[i] = min(effects)
		if(max(effects) > 0 & min(effects) < 0) {
			both.indicator[i] = 1
		}
	}
	new.genes = genes[both.indicator != 1]
	new.max.effects = max.effects[both.indicator != 1]
	new.min.effects = min.effects[both.indicator != 1]
	abs.min.effects = abs(new.min.effects)
	final.max.effects = vector('numeric', length(new.genes))
	for(i in 1:length(final.max.effects)) {
		if(new.max.effects[i] >= abs.min.effects[i]) {
			final.max.effects[i] = new.max.effects[i]
		} else if(absolute_effects == "TRUE") {
			final.max.effects[i] = abs.min.effects[i]
		} else {
			final.max.effects[i] = new.min.effects[i]
		}
	}
	
## Create genelist

	if(effect_type == "mean") {
		final.genes = genes
		final.effects = mean.effects
	} else if(effect_type == "max") {
		final.genes = new.genes
		final.effects = final.max.effects
	} else {
		stop("Invalid effect_type.")
	}
	final.effects[final.effects == 0] = min(final.effects[final.effects > 0])
	genelist = data.frame(cbind(final.genes, final.effects), stringsAsFactors = F)
	genelist[[2]] = as.numeric(genelist[[2]])
	colnames(genelist) = c('Gene', 'Effect')
	genelist = genelist[order(genelist[[2]], decreasing = T),]
	dataset.genes = as.character(genelist[,1])
	dataset.effects = genelist[,2]
	dataset.genes.N = length(dataset.genes)
	dataset.indexes = c(1:dataset.genes.N)
	
## Curate geneset database and remove genesets smaller than geneset_size_threshold

	print("Curating geneset database")

	geneset.matrix = as.matrix(read.table(geneset_database, row.names = 1, sep = '\t'))
	geneset.N = nrow(geneset.matrix)
	geneset.sizes = rowSums(geneset.matrix != 'null')
	max.geneset.size = max(geneset.sizes)
	geneset.names = rownames(geneset.matrix)

	hit.geneset.matrix = geneset.matrix
	miss.geneset.matrix = geneset.matrix
	geneset.hit.counts = vector(length = nrow(hit.geneset.matrix), mode = 'numeric')
	geneset.miss.counts = vector(length = nrow(miss.geneset.matrix), mode = 'numeric')
	geneset.hit.vector = vector(length = nrow(hit.geneset.matrix), mode = 'character')
	for(i in 1:nrow(geneset.matrix)) {
		hit.indicator = sign(match(geneset.matrix[i,], dataset.genes, nomatch = 0))
		hit.genes = geneset.matrix[i,][which(hit.indicator == 1)]
		miss.genes = geneset.matrix[i,][which(hit.indicator == 0)]
		hit.geneset.matrix[i,] = c(hit.genes, rep('null', times = ncol(geneset.matrix) - length(hit.genes)))
		miss.geneset.matrix[i,] = c(miss.genes, rep('null', times = ncol(geneset.matrix) - length(miss.genes)))
		geneset.hit.counts[i] = length(hit.genes)
		geneset.miss.counts[i] = length(miss.genes)
		geneset.hit.vector[i] = paste(hit.genes, collapse = ";")
	}

	geneset.size.indicator = which(geneset.hit.counts > geneset_size_threshold)
	geneset.matrix.2 = hit.geneset.matrix[geneset.size.indicator,]
	geneset.hit.vector.2 = geneset.hit.vector[geneset.size.indicator]
	geneset.names.2 = geneset.names[geneset.size.indicator]
	geneset.sizes.2 = geneset.sizes[geneset.size.indicator]
	geneset.N.2 = length(geneset.size.indicator)
	
## Compute enrichment scores for each geneset in original and permuted gene lists

	print("Computing enrichment scores")
	
	obs.ES = vector(length = geneset.N.2, mode = "numeric")
	obs.arg.ES = vector(length = geneset.N.2, mode = "numeric")
	obs.RES = matrix(nrow = geneset.N.2, ncol = dataset.genes.N)
	obs.indicator = matrix(nrow = geneset.N.2, ncol = dataset.genes.N)
	perm.ES = matrix(nrow = geneset.N.2, ncol = nperm)

	for(i in 1:geneset.N.2) {
		geneset = geneset.matrix.2[i,geneset.matrix.2[i,] != "null"]
		geneset.indexes = vector(length = length(geneset), mode = "numeric")
		geneset.indexes = match(geneset, dataset.genes)
		GSEA.results = GSEA_ES(gene.list = dataset.indexes, gene.set = geneset.indexes, weight = weight, correl.vector = dataset.effects, absolute_effects = absolute_effects)
		obs.ES[i] = signif(GSEA.results$ES, 3)
		obs.arg.ES[i] = GSEA.results$arg.ES
		obs.RES[i,] = GSEA.results$RES
		obs.indicator[i,] = GSEA.results$indicator
		for(r in 1:nperm) {
			new.dataset.indexes = sample(1:dataset.genes.N)
			GSEA.results = GSEA_ES(gene.list = new.dataset.indexes, gene.set = geneset.indexes, weight = weight, correl.vector = dataset.effects, absolute_effects = absolute_effects)
			perm.ES[i,r] = signif(GSEA.results$ES, 3)
		}
	}
	
## Normalize observed and permuted enrichment scores

	print("Normalizing enrichment scores")
	
	perm.ES.norm = matrix(nrow = geneset.N.2, ncol = nperm)
	obs.ES.norm = vector(length = geneset.N.2, mode = "numeric")
	
	for(i in 1:geneset.N.2) {
		mean.perm.ES = mean(abs(as.numeric(perm.ES[i,])))
		perm.ES.norm[i,] = signif(perm.ES[i,] / mean.perm.ES, 3)
		obs.ES.norm[i] = signif(obs.ES[i] / mean.perm.ES, 3)
	}

## Compute FDRs 

	print("Computing FDR q-values")

	FDR.ES = vector(length=geneset.N.2, mode="numeric")

	for(i in 1:geneset.N.2) {
			FDR.ES[i] = signif(sum(abs(perm.max.ES.norm[i,]) >= abs(obs.max.ES.norm[i])) / sum(abs(perm.max.ES.norm[i,]) >= 0), 3) # compute the fraction of ES scores larger than that observed for the original gene list across all permutations
	}

## Produce report and running enrichment plot for each gene set passing the FDR q-value cut-off

	print("Producing result tables and plots")

	time <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S") 
	for(i in 1:geneset.N.2) {
		
## Produce gene report
		
		gene.number = vector(length = geneset.sizes.2[i], mode = "character")
		gene.name = vector(length = geneset.sizes.2[i], mode = "character")
		gene.list.loc = vector(length = geneset.sizes.2[i], mode = "numeric")
		gene.effect = vector(length = geneset.sizes.2[i], mode = "numeric")
		gene.RES = vector(length = geneset.sizes.2[i], mode = "numeric")
		hairpin.index = vector(length = geneset.sizes.2[i], mode = "character")
		hairpin.count = vector(length = geneset.sizes.2[i], mode = "numeric")
		hairpin.effects = vector(length = geneset.sizes.2[i], mode = "character")
		hairpin.sd = vector(length = geneset.sizes.2[i], mode = "numeric")

		kk = 1
		if(obs.ES[i] >= 0) {
			order = seq(1, dataset.genes.N, 1)
		} else {
			order = seq(dataset.genes.N, 1, -1)
		}
		for(k in order) {
			if(obs.indicator[i, k] == 1) {
				gene.number[kk] = kk
				gene.name[kk] = dataset.genes[k]
				gene.list.loc[kk] = k
				gene.effect[kk] = dataset.effects[k]
				hairpin.index[kk] = paste(which(data[[2]] == dataset.genes[k]), collapse = ';')
				hairpin.count[kk] = length(strsplit(hairpin.index[kk], ';')[[1]])
				hairpin.effects[kk] = paste(data[as.numeric(strsplit(hairpin.index[kk], ';')[[1]]),3], collapse = ';')
				hairpin.sd[kk] = signif(sd(as.numeric(strsplit(hairpin.effects[kk], ';')[[1]])), 3)
				kk = kk + 1
			}
		}
		
		gene.report = data.frame(cbind(gene.number, gene.name, gene.list.loc, gene.effect, hairpin.count, hairpin.effects, hairpin.sd))
		names(gene.report) = c("#", "GENE", "LIST LOCATION", "EFFECT", 'HAIRPINS', 'HAIRPIN EFFECTS', 'HAIRPIN SD')
		gene.report = gene.report[gene.report[[2]] != '',]
		filename = paste(output.directory, time, '_', substring(geneset.names.2[i], 1, 40), ".txt", sep = "", collapse = "")
		write.table(gene.report, file = filename, quote = F, row.names = F, sep = "\t")
		
## Produce running enrichment plot
		
		pdf.filename = paste(output.directory, time, '_', substring(geneset.names.2[i], 1, 40), ".pdf", sep = "", collapse = "")
		pdf(file = pdf.filename, height = 6, width = 14)
		ind = 1:dataset.genes.N
		min.RES = min(obs.RES[i,])
		max.RES = max(obs.RES[i,])
		if(max.RES < 0.3) max.RES = 0.3
		if(min.RES > -0.3) min.RES = -0.3
		delta = (max.RES - min.RES) * 0.50
		min.plot = min.RES - 2 * delta
		max.plot = max.RES
		max.corr = max(dataset.effects)
		min.corr = min(dataset.effects)
		obs.correl.vector.norm = (dataset.effects - min.corr) / (max.corr - min.corr) * 1.25 * delta + min.plot
		zero.corr.line = (-min.corr / (max.corr - min.corr)) * 1.25 * delta + min.plot
		sub.string = paste("Number of genes: ", dataset.genes.N, " (in list), ", geneset.sizes.2[i], " (in gene set)", sep = "", collapse = "")
		main.string = geneset.names.2[i]
		plot(ind, obs.RES[i,], main = main.string, sub = sub.string, xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)", xlim = c(1, dataset.genes.N), ylim = c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, cex.main = 1, col = 2)
		for(j in seq(1, dataset.genes.N, 10)) {
			lines(c(j, j), c(zero.corr.line, obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
		}
		lines(c(1, dataset.genes.N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
		lines(c(obs.arg.ES[i], obs.arg.ES[i]), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 2) # max enrichment vertical line
		for(j in 1:dataset.genes.N) {
			if(obs.indicator[i, j] == 1) {
				lines(c(j, j), c(min.plot + 1.25 * delta, min.plot + 1.75 * delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
			}
		}
		lines(ind, obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
		lines(c(1, dataset.genes.N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
		if(ES_type == 'max') {
			leg.txt = paste("Max at ", obs.arg.ES[i], sep = "", collapse = "")
		} else {
			leg.txt = paste("Med at ", obs.arg.ES[i], sep = "", collapse = "")
		}
		text(x = obs.arg.ES[i], y = min.plot + 2 * delta, labels = leg.txt, cex = 1.0)
		dev.off()
	}
	
## Produce global results report

	hit.count = rowSums(obs.indicator)
	hit.percent = hit.count / geneset.sizes.2
	global.report = data.frame(cbind(geneset.names.2, geneset.sizes.2, geneset.hit.vector.2, obs.ES, obs.ES.norm, FDR.ES),stringsAsFactors=F)
	names(global.report) = c("GENESET NAME", "SIZE", 'GENES', "ES", "NES", "FDR")
	global.report = global.report[order(global.report$"FDR"),]
	class(global.report$"FDR") = 'numeric'
	header = as.vector(c(
		paste('input = ', input, sep = ''),
		paste('annotation_file = ', annotation_file, sep = ''),
		paste('geneset_database = ', geneset_database, sep = ''),
		paste('pval_threshold = ', pval_threshold, sep = ''),
		paste('fdr_threshold = ', fdr_threshold, sep = ''),
		paste('absolute_effects = ', absolute_effects, sep = ''),
		paste('effect_type = ', effect_type, sep = ''),
		paste('ES_type = ', ES_type, sep = ''),
		paste('geneset_size_threshold = ', geneset_size_threshold, sep = ''),
		paste('nperm = ',nperm,sep=''),
		paste('weight = ',weight,sep=''),''),
	mode = 'character')
	global.filename = paste(output.directory, time, '_', "global_report.txt", sep = "", collapse = "")
	writeLines(header, con = global.filename, sep = '\n')
	write.table(global.report, file = global.filename, quote = F, row.names = F, sep = "\t", append = T)
	
## Produce summary histograms

	hits = global.report[global.report$"FDR" <= fdr_threshold, 1]
	total = global.report[global.report$"FDR" <= 1, 1]
	hit.matrix = matrix(nrow = length(hits), ncol = ncol(geneset.matrix.2))
	total.matrix = matrix(nrow = length(total), ncol = ncol(geneset.matrix.2))
	for(i in 1:length(hits)) {
		hit.matrix[i,] = geneset.matrix.2[which(geneset.names.2 == hits[i]),]
	}
	for(i in 1:length(total)) {
		total.matrix[i,] = geneset.matrix.2[which(geneset.names.2 == total[i]),]
	}
	hit.vector = hit.matrix[1,]
	total.vector = total.matrix[1,]
	for(i in 2:nrow(hit.matrix)) {
		hit.vector = c(hit.vector,hit.matrix[i,])
	}
	for(i in 2:nrow(total.matrix)) {
		total.vector = c(total.vector,total.matrix[i,])
	}
	hit.vector = hit.vector[hit.vector != 'null']
	hit.table = table(hit.vector)
	hit.table = hit.table[order(hit.table, decreasing = T)]
	total.vector = total.vector[total.vector != 'null']
	total.table = data.frame(table(total.vector))
	denom.table = hit.table
	effect.table = hit.table
	for(i in 1:length(hit.table)) {
		denom.table[i] = total.table[which(total.table[,1] == names(hit.table[i])),2]
		effect.table[i] = genelist[which(genelist[,1] == names(hit.table[i])),2]
	}
	hit.table = hit.table[order(effect.table, decreasing = T)]
	denom.table = denom.table[order(effect.table, decreasing = T)]
	effect.table = effect.table[order(effect.table, decreasing = T)]
	ratio.table = hit.table / denom.table
	barplot.filename = paste(output.directory, time, '_', "histogram", ".pdf", sep = "", collapse = "")
	pdf(file = barplot.filename, height = 14, width = 30)
	par(mfcol = c(2,1), mar = c(5,4,4,2))
	xpos = barplot(ratio.table[hit.table > 1], space = .5, cex.names = .8, las = 2, main = paste('Fraction of Genesets with FDR <= ', fdr_threshold, sep = ''))
	text(x = xpos, y = ratio.table[hit.table > 1] + .1 * max(ratio.table[hit.table > 1]), labels = sub('0.', '.', round(ratio.table[hit.table > 1], 1)), xpd = TRUE, cex = .7)
	barplot(effect.table[hit.table > 1], space = .5, cex.names = .8, las = 2, main = 'Effect Size')
	text(x = xpos, y = effect.table[hit.table > 1] + .1 * max(effect.table[hit.table > 1]), labels = sub('0.', '.', round(effect.table[hit.table > 1], 1)), xpd = TRUE, cex = .7)
	dev.off()
}

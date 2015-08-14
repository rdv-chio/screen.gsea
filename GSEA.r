GSEA = function(input, annotation_file, geneset_database, output.directory, pval_threshold = 1, fdr_threshold = .1, absolute_effects = "TRUE", effect_type = "med", ES_type = "med", geneset_size_threshold = 3, nperm = 1000, weight = 1) {

## input: Input filename. Columns of input should be: shrna, gene, effect, pval.
## annotation_file: Current gene annotations for each hairpin. Columns should be: shrna, gene.
## geneset_database: Genesets to be analyzed.
## output.directory: Directory where results should be stored.
## pval_threshold: Numeric 0-1. Hairpins with p-values greater than this threshold will be removed from the dataset.
## fdr_threshold: Numeric 0-1. Genesets with FDR values greater than this threshold will be excluded from summary plots.
## absolute_effects: Logical. If "TRUE", all effects will be converted to positive values. 
## effect_type: "med" or "max". Determines whether mean or median effect should be used for each gene.
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
	
## Calculate median and max effect for each gene

	genes = unique(data[[2]])
	med.effects = vector('numeric', length(genes))
	max.effects = vector('numeric', length(genes))
	min.effects = vector('numeric', length(genes))
	final.effects = vector('numeric', length(genes))
	
	for(i in 1:length(genes)) {
		effects = data[data[[2]] == genes[i],3]
		med.effects[i] = median(effects)
		max.effects[i] = max(effects)
		min.effects[i] = min(effects)
	}
	abs.min.effects = abs(min.effects)
	abs.med.effects = abs(med.effects)
	
## Create genelist
	
	if(absolute_effects == "TRUE" & effect_type == 'med') {
		final.effects = abs.med.effects
	} else if(absolute_effects == "FALSE" & effect_type == 'med') {
		final.effects = med.effects
	} else if(effect_type == 'max') {
		for(i in 1:length(final.effects)) {
			if(max.effects[i] >= abs.min.effects[i]) {
				final.effects[i] = max.effects[i]
			} else if(absolute_effects == "TRUE") {
				final.effects[i] = abs.min.effects[i]
			} else {
				final.effects[i] = min.effects[i]
			}
		}
	} else {
		stop("Invalid effect_type.")
	}

	final.effects[final.effects == 0] = min(final.effects[final.effects > 0])
	genelist = data.frame(cbind(genes, final.effects), stringsAsFactors = F)
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
	
	perm.matrix = replicate(nperm, sample(1:dataset.genes.N))
	perm.ES = matrix(NA, nrow = geneset.N.2, ncol = nperm)
	correl.vector = abs(dataset.effects ** weight)
	geneset.rows = lapply(split(geneset.matrix.2, row(geneset.matrix.2)), function(x) x[x != 'null'])
	geneset.index.list = lapply(geneset.rows, match, table = dataset.genes)
	obs.GSEA.results = lapply(geneset.index.list, GSEA_ES, gene.list = dataset.indexes, correl.vector = correl.vector, absolute_effects = absolute_effects, ES_type = ES_type)
	obs.ES = sapply(obs.GSEA.results, "[[", 1)
	obs.arg.ES = sapply(obs.GSEA.results, "[[", 2)
	obs.RES = matrix(sapply(obs.GSEA.results, "[[", 3), nrow = geneset.N.2, ncol = length(dataset.indexes), byrow=T)
	obs.indicator = matrix(sapply(obs.GSEA.results, "[[", 4), nrow = geneset.N.2, ncol = length(dataset.indexes), byrow=T)
	if(ES_type == "max") {	
		for(i in 1:geneset.N.2) {
			perm.ES[i,] = sapply(split(perm.matrix, col(perm.matrix)), GSEA_ES_perm_max, geneset.index.list[[i]], correl.vector = correl.vector)
		}
	} else if(ES_type == 'med') {
		for(i in 1:geneset.N.2) {
			perm.ES[i,] = sapply(split(perm.matrix, col(perm.matrix)), GSEA_ES_perm_med, geneset.index.list[[i]], correl.vector = correl.vector)
		}
	} else {
		stop("Invalid ES_type.")
	} 
	
## Normalize observed and permuted enrichment scores

	print("Normalizing enrichment scores")
	
	perm.ES.norm = matrix(nrow = geneset.N.2, ncol = nperm)
	obs.ES.norm = vector(length = geneset.N.2, mode = "numeric")
	
	for(i in 1:geneset.N.2) {
		mean.perm.ES = mean(abs(as.numeric(perm.ES[i,])))
		perm.ES.norm[i,] = perm.ES[i,] / mean.perm.ES
		obs.ES.norm[i] = obs.ES[i] / mean.perm.ES
	}

## Compute FDRs 

	print("Computing FDR q-values")

	FDR.ES = vector(length=geneset.N.2, mode="numeric")

	for(i in 1:geneset.N.2) {
			FDR.ES[i] = sum(abs(perm.ES.norm[i,]) >= abs(obs.ES.norm[i])) / nperm # compute the fraction of ES scores larger than that observed for the original gene list across all permutations
	}

## Produce geneset report and running enrichment plot for each gene set passing the FDR q-value cut-off

	print("Producing result tables and plots")

	time = format(Sys.time(), "%Y_%m_%d_%H_%M_%S") 
	
	for(i in 1:geneset.N.2) {	
	
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
			if(obs.indicator[i,k] == 1) {
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
		
		pdf.filename = paste(output.directory, time, '_', substring(geneset.names.2[i], 1, 40), ".pdf", sep = "", collapse = "")
		pdf(file = pdf.filename, height = 6, width = 6)
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
			if(obs.indicator[i,j] == 1) {
				lines(c(j,j), c(min.plot + 1.25 * delta, min.plot + 1.75 * delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
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
	
## Write parameters report

	parameters = as.vector(c(
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
	parameters.filename = paste(output.directory, time, '_', "parameters.txt", sep = "", collapse = "")
	writeLines(parameters, con = parameters.filename, sep = '\n')
	
## Write global results report

	obs.ES = signif(obs.ES, 3)
	obs.ES.norm = signif(obs.ES.norm, 3)
	FDR.ES = signif(FDR.ES, 3)
	
	global.report = data.frame(cbind(geneset.names.2, geneset.sizes.2, geneset.hit.vector.2, obs.ES, obs.ES.norm, FDR.ES),stringsAsFactors=F)
	names(global.report) = c("GENESET NAME", "SIZE", 'GENES', "ES", "NES", "FDR")
	global.report = global.report[order(global.report$"FDR"),]
	class(global.report$"FDR") = 'numeric'
	global.filename = paste(output.directory, time, '_', "global_report.txt", sep = "", collapse = "")
	write.table(global.report, file = global.filename, quote = F, row.names = F, sep = "\t", append = T)
	
## Draw summary histograms

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
	total.vector = total.vector[total.vector != 'null']
	hit.table = table(hit.vector)
	total.table = table(total.vector)
	denom.table = hit.table
	effect.table = hit.table
	for(i in 1:length(hit.table)) {
		denom.table[i] = total.table[which(names(total.table) == names(hit.table[i]))]
		effect.table[i] = genelist[which(genelist[,1] == names(hit.table[i])),2]
	}
	p.table = denom.table / geneset.N.2
	invp.table = 1 - p.table
	pk.table = p.table ** hit.table
	nk.table = length(hits) - hit.table
	binomial.table = pk.table * invp.table ** nk.table
	p.threshold = .05 / length(binomial.table)
	delta = p.threshold - min(binomial.table)
	order = order(binomial.table)
	binomial.table = binomial.table[order]
	effect.table = effect.table[order]
	subset = which(binomial.table <= p.threshold + delta)
	binomial.table = binomial.table[subset]
	effect.table = effect.table[subset]
	barplot.filename = paste(output.directory, time, '_', "histogram", ".pdf", sep = "", collapse = "")
	pdf(file = barplot.filename, height = 14, width = length(binomial.table)/8)
	par(mfcol = c(2,1), mar = c(5,4,4,2))
	xpos = barplot(binomial.table, space = .5, cex.names = .8, axes = F, las = 2, ylim = c(min(binomial.table), p.threshold + delta), main = paste('Binomial Probability of N Genesets with FDR <= ', fdr_threshold, sep = ''))
	abline(h = p.threshold)
	axis(side = 1, line = 0, tick = F, labels = F)
	axis(side = 2, line = 0)
	barplot(effect.table, space = .5, cex.names = .8, axes = F, las = 2, main = 'Effect Size')
	text(x = xpos, y = effect.table + .025 * max(effect.table), labels = sub('0.', '.', round(effect.table, 1)), xpd = TRUE, cex = .7)
	axis(side = 1, line = 0, tick = F, labels = F)
	axis(side = 2, line = 0)
	dev.off()
	
## Draw summary heatmap

	library(pheatmap)
	heatmap.genes = names(binomial.table[binomial.table <= p.threshold])
	heatmap.effects = effect.table[which(binomial.table <= p.threshold)]
	hit.fdrs = global.report[global.report$"FDR" <= fdr_threshold, 6]
	hit.fdrs[hit.fdrs == 0] = min(hit.fdrs[hit.fdrs > 0])
	heatmap.fdrs = -log(hit.fdrs, 10)
	corr.matrix = matrix(NA,nrow = nrow(hit.matrix), ncol = length(heatmap.genes))
	for(i in 1:nrow(corr.matrix)) {
			corr.matrix[i,] = sign(match(heatmap.genes,hit.matrix[i,],nomatch=0))
	}
	rownames(corr.matrix) = substring(hits, 1, 20)
	colnames(corr.matrix) = heatmap.genes
	for(i in 1:nrow(corr.matrix)) {
		for(j in 1:ncol(corr.matrix)) {
			corr.matrix[i,j] = corr.matrix[i,j] * heatmap.effects[j] * heatmap.fdrs[i]
		}
	}
	corr.matrix[corr.matrix == 0] = NA
	corr.matrix[!is.na(corr.matrix)] = log(corr.matrix[!is.na(corr.matrix)], 10)
	corr.matrix[is.na(corr.matrix)] = min(corr.matrix[!is.na(corr.matrix)]) - 1
	corr.matrix = corr.matrix - min(corr.matrix)
	heatmap.filename = paste(output.directory, time, '_', 'heatmap', ".pdf", sep = "", collapse = "")
	lower_limit = 0.9 * min(corr.matrix[corr.matrix > 0])
	upper_limit = 1.1 * max(corr.matrix)
	delta = upper_limit - lower_limit
	pheatmap(corr.matrix, color = c('grey100', hsv(h = 0, s = seq(0,1,1/20), v = seq(1,.8,-.2/20))), breaks = c(0, seq(lower_limit, upper_limit, delta / 19)), border_color = 'grey50', treeheight_row = 0, treeheight_col = 0, legend = T, filename = heatmap.filename, height = length(hits) / 5, width = length(heatmap.genes) / 7)
}

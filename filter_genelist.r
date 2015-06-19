## to filter a gene list based on pval and assign one effect value per gene

filter_genelist = function(input,annotation_file,pval_threshold,absolute_effects="FALSE",effect_type="max",output) {

## input: Input filename. Columns of input should be: shrna, gene, effect, pval.
## annotation_file: File containing current gene annotations for each hairpin. Columns should be: shrna, gene.
## pval_threshold: Numeric 0-1. Hairpins with p-values greater than this threshold will be removed from the dataset.
## absolute_effects: Logical. If "TRUE", all effects will be converted to positive values. 
## effect_type: "mean" or "max". Determines whether mean or maximum effect should be used for each gene.
## output: Output filename.

## read dataset

	data = read.table(input,header=T,sep='\t',colClasses=c('character','character','numeric','numeric'))
	
## update shRNA annotations
	
	annotation_file = read.table(annotation_file,header=T,sep='\t',colClasses=c('character','character'))
	annotations = vector('character',length=length(data[[1]]))
	for(i in 1:length(data[[1]])) {
		annotations[i] = annotation_file[annotation_file[[1]] == data[[1]][i],2]
	}
	data = data.frame(cbind(data,annotations),stringsAsFactors=F)
	data = data[,c(1,5,3,4)]
	data[[3]] = as.numeric(data[[3]])
	data[[4]] = as.numeric(data[[4]])
	colnames(data) = c('shRNA','Gene','Effect','P Value')

## filter based on pval_threshold
	
	data = data[data[[4]]<=pval_threshold,]
	
## calculate mean effect for each gene

	genes = unique(data[[2]])
	mean.effects = vector('numeric',length(genes))
	for(i in 1:length(genes)) {
		effects = data[data[[2]] == genes[i],3]
		mean.effects[i] = mean(effects)
	}
	if(absolute == "TRUE") {
		mean.effects = abs(mean.effects)
	} else {
	}
		
## calculate max effect for each gene (genes with both positive and negative effects are removed)

	both.indicator = vector('numeric',length(genes))
	max.effects = vector('numeric',length(genes))
	min.effects = vector('numeric',length(genes))
	for(i in 1:length(genes)) {
		effects = data[data[[2]] == genes[i],3]
		max.effects[i] = max(effects)
		min.effects[i] = min(effects)
		if(max(effects)>0 & min(effects)<0) {
			both.indicator[i] = 1
		}
	}
	new.genes = genes[both.indicator != 1]
	new.max.effects = max.effects[both.indicator != 1]
	new.min.effects = min.effects[both.indicator != 1]
	abs.min.effects = abs(new.min.effects)
	final.max.effects = vector('numeric',length(new.genes))
	for(i in 1:length(final.max.effects)) {
		if(new.max.effects[i]>=abs.min.effects[i]) {
			final.max.effects[i] = new.max.effects[i]
		} else if(absolute == "TRUE") {
			final.max.effects[i] = abs.min.effects[i]
		} else {
			final.max.effects[i] = new.min.effects[i]
		}
	}
	
## write output

	if(effect_type == "mean") {
		final.genes = genes
		final.effects = mean.effects
	} else if(effect_type == "max") {
		final.genes = new.genes
		final.effects = final.max.effects
	} else {
		stop("Invalid effect_type.")
	}
	final.effects[final.effects==0] = min(final.effects)
	new.data = data.frame(cbind(final.genes,final.effects),stringsAsFactors=F)
	new.data[[2]] = as.numeric(new.data[[2]])
	colnames(new.data) = c('Gene','Effect')
	new.data = new.data[order(new.data[[2]],decreasing=T),]
	write.table(new.data,output,row.names=F,sep='\t')
}

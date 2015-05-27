## to filter a gene list based on pval and assign one (max effect) hairpin per gene
## columns of input should be: shrna, gene, effect, pval

filter_genelist = function(input,pval_threshold,output) {

## read dataset and filter based on pval_threshold

	data = read.table(input,header=T,sep='\t',colClasses=c('character','factor','numeric','numeric'))
	data = data[data[[4]]<pval_threshold,]
	
## remove genes with hairpins with both positive and negative effects

	genes = levels(data[[2]])
	both.indicator = vector('numeric',length(genes))
	max.effects = vector('numeric',length(genes))
	min.effects = vector('numeric',length(genes))
	for(i in 1:length(both.indicator)) {
		effects = data[data[[2]] == genes[i],3]
		max.effects[i] = max(effects)
		min.effects[i] = min(effects)
		if(max(effects)>0 & min(effects)<0) {
			both.indicator[i] = 1
		}
	}
	new.genes = genes[both.indicator != 1]
	
## find max effect for each gene
	
	new.max.effects = max.effects[both.indicator != 1]
	new.min.effects = min.effects[both.indicator != 1]
	abs.min.effects = abs(new.min.effects)
	final.effects = vector('numeric',length(new.genes))
	for(i in 1:length(final.effects)) {
		if(new.max.effects[i]>=abs.min.effects[i]) {
			final.effects[i] = new.max.effects[i]
		} else {
			final.effects[i] = new.min.effects[i]
		}
	}
	new.data = data.frame(cbind(new.genes,max.effects),stringsAsFactors=F)
	new.data[[2]] = as.numeric(new.data[[2]])
	colnames(new.data) = c('Gene','Max Effect')
	new.data = new.data[order(new.data[[2]],decreasing=T),]
	write.table(new.data,output,row.names=F,sep='\t')
}

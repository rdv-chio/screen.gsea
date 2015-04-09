## to find the hairpin with the largest effect size for each gene given a screen dataset with mpd

find_max_hairpin = function(input,output) {
	data = read.table(input,header=T,sep='\t',colClasses=c('character','factor','numeric','numeric'))
	genes = levels(data[[2]])
	both.indicator = vector('numeric',length(genes))
	for(i in 1:length(both.indicator)) {
		effects = data[data[[2]] == genes[i],3]
		if(max(effects)>0 & min(effects)<0) {
			both.indicator[i] = 1
		}
	}
	new.genes = genes[both.indicator != 1]
	max.effects = vector('numeric',length(new.genes))
	for(i in 1:length(max.effects)) {
		if(mean(data[data[[2]] == new.genes[i],3])>0) {
			max.effects[i] = max(data[data[[2]] == new.genes[i],3])
		} else {
			max.effects[i] = min(data[data[[2]] == new.genes[i],3])
		}
	}
	new.data = data.frame(cbind(new.genes,max.effects),stringsAsFactors=F)
	new.data[[2]] = as.numeric(new.data[[2]])
	colnames(new.data) = c('Gene','Max Effect')
	new.data = new.data[order(new.data[[2]],decreasing=T),]
	write.table(new.data,output,row.names=F,sep='\t')
	new.data
}
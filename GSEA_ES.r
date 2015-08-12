GSEA_ES <- function(gene.list, gene.set, weight = 1, correl.vector = NULL, absolute_effects = "TRUE") {  

# Inputs:
#   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
#   weight: Exponential scalar for effect sizes
#   correl.vector: A vector with the effect sizes corresponding to the genes in the gene list 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
#   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
#   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 

	tag.indicator = sign(match(gene.list, gene.set, nomatch=0))
	no.tag.indicator = 1 - tag.indicator 
	Ngl = length(gene.list)
	Ngs = length(gene.set) 
	Ndiff =  Ngl - Ngs 
	correl.vector = abs(correl.vector ** weight)
	norm.tag = 1.0 / sum(correl.vector[tag.indicator == 1]) # normalizing so that genesets with lots of hits dont have artificially high ES scores; sum(correl.vector*norm.tag)=1
	norm.no.tag = 1.0 / Ndiff # if geneset is small compared to gene list, decay is slow. if geneset is large, decay is faster
	RES = cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
	max.ES = signif(max(RES), 3)
	arg.max.ES = which.max(RES)
	min.ES = signif(min(RES), 3)
	arg.min.ES = which.min(RES)
	med.ES = signif(median(RES[tag.indicator == 1]), 3)
	arg.med.ES = median(which(tag.indicator == 1))
	
	if(absolute_effects == "TRUE") {
		if(ES_type == "max") {
			ES = max.ES
			arg.ES = arg.max.ES
		} else {
			ES = med.ES
			arg.ES = arg.med.ES
		}
	} else {
		if(ES_type == "max") {
			if (max.ES > -min.ES) {
				ES = max.ES
				arg.ES = arg.max.ES
			} else {
				ES = min.ES
				arg.ES = arg.min.ES
			} 
		} else {
			ES = med.ES
			arg.ES = arg.med.ES
		}
	}
	return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}
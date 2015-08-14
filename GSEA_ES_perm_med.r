GSEA_ES_perm_med <- function(gene.list, gene.set, correl.vector) {  

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
	norm.tag = 1.0 / sum(correl.vector[tag.indicator == 1]) # normalizing so that genesets with lots of hits dont have artificially high ES scores; sum(correl.vector*norm.tag)=1
	norm.no.tag = 1.0 / Ndiff # if geneset is small compared to gene list, decay is slow. if geneset is large, decay is faster
	RES = cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
	ES = median(RES[tag.indicator == 1])
	return(ES)
}
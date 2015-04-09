GSEA_ES_lite_2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  

# Inputs:
#   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
#   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
#   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
#   correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
#
# Outputs:
#   ES: Enrichment score (real number between -1 and +1) 
#   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
#   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
#   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 

	tag.indicator <- sign(match(gene.list, gene.set, nomatch=0)) # 0 (no tag) or 1 (tag) 
	no.tag.indicator <- 1 - tag.indicator 
	N <- length(gene.list) 
	Nh <- length(gene.set) 
	Nm <-  N - Nh 
	if (weighted.score.type == 0) {
		correl.vector <- rep(1, N)
	}
	alpha <- weighted.score.type
	correl.vector <- abs(correl.vector**alpha)
	sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
	norm.tag <- 1.0/sum.correl.tag # normalizing so that sum(correl.vector * norm.tag) is closer to 1.0
	norm.no.tag <- 1.0/Nm # if geneset is small compared to gene list, decay is slow. if geneset is large, decay is faster
	RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
	max.ES <- max(RES)
	min.ES <- min(RES)
	if (max.ES > - min.ES) {
		ES <- signif(max.ES, digits = 5)
		arg.ES <- which.max(RES)
	} else {
		ES <- signif(min.ES, digits=5)
		arg.ES <- which.min(RES)
	}
	return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}
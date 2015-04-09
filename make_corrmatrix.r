## to make geneset correlation matrix given a geneset database (gmt file)

make_corrmatrix = function(input,output) {
	genesets <- readLines(input)
	geneset.N <- length(genesets)
	geneset.sizes <- vector(length = geneset.N, mode = "numeric") 
	for (i in 1:geneset.N) {
		geneset.sizes[i] <- length(unlist(strsplit(genesets[[i]], "\t"))) - 2
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

	corr.matrix = matrix(NA,nrow = nrow(geneset.matrix), ncol = length(refseq))
	for(i in 1:nrow(geneset.matrix)) {
		for(j in 1:length(refseq)) {
			corr.matrix[i,j] = ifelse(is.element(refseq[j],geneset.matrix[i,]),1,0)
		}
	}
	rownames(corr.matrix) = geneset.names
	colnames(corr.matrix) = refseq
	corr.matrix = corr.matrix[,colSums(corr.matrix) != 0]
	write.table(corr.matrix,output,sep='\t')
}
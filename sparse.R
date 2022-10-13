#R functions to manipulate sparse matrix

#Convert a sparse matrix to a dense matrix.
#Use this rather than as.matrix when the matrix is too large that causes "Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 105"
#Ref: https://programmerah.com/the-sparse-matrix-of-r-language-is-too-large-to-be-used-as-matrix-8856/
as_matrix <- function(x, verbose = TRUE){
	t1 = Sys.time()
	if(inherits(x, 'dgCMatrix')){
		if(verbose) cat('Convert dgCMatrix matrix to dense matrix...\n')
		tmp <- matrix(data = 0, nrow = x@Dim[1], ncol = x@Dim[2])
		row_pos <- x@i + 1
		col_pos <- findInterval(seq(x@x) - 1, x@p[-1]) + 1
		val <- x@x
		for (i in seq_along(val)){
			tmp[row_pos[i], col_pos[i]] <- val[i]
		}
	}else if(inherits(x, 'dgTMatrix')){
		if(verbose) cat('Convert dgTMatrix matrix to dense matrix...\n')
		tmp <- matrix(data = 0, nrow = x@Dim[1], ncol = x@Dim[2])
		tmp[1 + x@i + (as.double(x@j) * x@Dim[1])] <- x@x
	}else{
		stop('Input x must be a dgCMatrix or dgTMatrix object\n')
	}
	rownames(tmp) <- x@Dimnames[[1]]
	colnames(tmp) <- x@Dimnames[[2]]
	if(verbose) cat('Conversion done after', round(Sys.time() - t1, 1), 'seconds.\n')
	return(tmp)
}

#Convert a sparse matrix to a dense matrix.
dgT2dgCMatrix <- function(x, verbose = TRUE){
	#Use this rather than as.matrix when the matrix is too large that causes "Cholmod error 'problem too large' at file ../Core/cholmod_dense.c, line 105"
	#https://programmerah.com/the-sparse-matrix-of-r-language-is-too-large-to-be-used-as-matrix-8856/
	t1 = Sys.time()
	stopifnot(inherits(x, 'dgTMatrix'))
	#
	#Function colNN counts the number of non-zero elements in each column
	#x is an integer vector with value 1 ~ dim
	#return a vector of counts, with first element of zero in the beginning
	#
	Rcpp:::cppFunction("
	IntegerVector colNN(IntegerVector x, int dim){
		int n = x.size();
		IntegerVector r(dim + 1);
		for (int i = 0; i < n; ++i) {
			r[x[i]] ++;
		}
		return r; 
	}
	")
	if(verbose) cat('Convert dgTMatrix to dgCMatrix...\n')
	tmp <- new('dgCMatrix')
	tmp@Dim <- x@Dim
	tmp@Dimnames <- x@Dimnames
	tmp@factors <- x@factors
	tmp@x <- numeric(length(x@x))
	tmp@i <- integer(length(x@x))
	tmp@p <- cumsum(colNN(x@j + 1, x@Dim[2]))
	# tmp@p <- integer(x@Dim[2] + 1L)
	# #for(j in x@j) tmp@p[j + 2] <- tmp@p[j + 2] + 1L #slow in R, use colNN
	# nJ <- table(sort(x@j))
	# nJn <- as.integer(names(nJ))
	# nJ <- as.vector(nJ)
	# tmp@p[1] <- 0L
	# tmp@p[nJn + 2] <- nJ
	# tmp@p <- cumsum(tmp@p)
	ord <- order(x@j, x@i)
	tmp@i <- x@i[ord]
	tmp@x <- x@x[ord]
	if(verbose) cat('Conversion done after', round(Sys.time() - t1, 1), 'seconds.\n')
	return(tmp)
}

#cbind two dgCmatrices
cbind_dgCMatrices <- function(x, y){
	stopifnot(inherits(x, 'dgCMatrix'))
	stopifnot(inherits(y, 'dgCMatrix'))
	stopifnot(x@Dim[1] == y@Dim[1])
	stopifnot(identical(x@Dimnames[[1]], y@Dimnames[[1]]))
	if(is.null(x@Dimnames[[2]])){
		if(! is.null(y@Dimnames[[2]])) stop('x has colnames but y has not')
	}else{
		if(is.null(y@Dimnames[[2]])) stop('y has colnames but x has not')
		x@Dimnames[[2]] <- c(x@Dimnames[[2]], y@Dimnames[[2]])
	}
	x@i <- c(x@i, y@i)
	x@p <- c(x@p, x@p[length(x@p)] + y@p[-1])
	x@x <- c(x@x, y@x)
	x@Dim[2] <- x@Dim[2] + y@Dim[2]
	x
}

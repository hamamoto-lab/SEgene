suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(openxlsx))
suppressMessages(library(pbmcapply))
suppressMessages(library(rhdf5))
suppressMessages(library(stringr))


load_data <- function(fnameE, fnameP, verbose=T) {
	#' Loads gene expression and peak calling files as dataframes
	#'
	#' @param fnameE character Gene expression file
	#' @param fnameP character Peak calling file
	#' @param verbose logical Verbose messages boolean
	#'
	#' @return List of filtered and ordered dataframes  

	#Loading the gene expression file as a dataframes
	if (verbose) {
		cat("loading", fnameE, "\n")
	}
	datE <- fread(fnameE, sep=",", stringsAsFactors=F, check.names=F)
	colnames(datE)[1:5] <- c("symbol", "chr", "strand", "start", "end")

	#Loading the peak calling file as a dataframes
	if (verbose) {
		cat("loading", fnameP, "\n")
	}
	datP <- fread(fnameP, stringsAsFactors=F, check.names=F)
	colnames(datP)[1:4] <- c("PeakID", "Chr", "Start", "End")
	
	#Getting gene expression and peak calling sample names
	sampleNameE <- colnames(datE)[-(1:5)]
	sampleNameP <- colnames(datP)[-(1:4)]

	#Getting sample names common in gene expression and peak calling datasets
	sampleName <- intersect(sampleNameE, sampleNameP)

	#Slice dataframes to only contain columns in common
	datE <- cbind(datE[, 1:5], datE[, ..sampleName], stringsAsFactors=F)
	datP <- cbind(datP[, 1:4], datP[, ..sampleName], stringsAsFactors=F)

	#Checking for duplicate genes in gene expression dataframe
	length(unique(datE$symbol)) == nrow(datE)

	#Remove genes on mitochondria, unknown chromosomes and sequences with unknown orientation
	chrs <- paste("chr", c(1:22, "X", "Y"), sep="")
	datE <- datE[is.element(datE$chr, chrs), ]
	datP <- datP[is.element(datP$Chr, chrs), ]

	#Ordering the gene expression dataframe
	datE <- datE[order(datE$start, datE$end), ]
	datE <- datE[str_order(datE$chr, numeric=T), ]

	#Ordering the peak calling dataframe
	datP <- datP[order(datP$Start, datP$End), ]
	datP <- datP[str_order(datP$Chr, numeric=T), ]

	return(list(datE=datE, datP=datP))
}


find_link <- function(dat, txWidth, min.var.quantile.E, min.var.quantile.P, seed, verbose) {
	#' Finds links between genes and peaks
	#'
	#' @param dat list List of filtered and ordered gene expression and peak calling dataframes
	#' @param txWidth integer Gene start search space range
	#' @param min.var.quantile.E double Bottom percentage of genes to remove from analysis
	#' @param min.var.quantile.P double Bottom percentage of peaks to remove from analysis
	#' @param seed integer Random seed number
	#' @param verbose logical Verbose messages boolean
	#'
	#' @return List of gene search space and peak links dataframes and Genomic Ranges objects

	#Initialising variables
	chrs <- paste("chr", c(1:22, "X", "Y"), sep="")

	if (verbose) {
		cat("finding links between peaks and transcriptional regulatory regions", "\n")
	}

	#Split peak calling dataframe into peak information and peak values dataframes
	datPa <- dat$datP[, c("Chr", "Start", "End", "PeakID")]
	datPd <- dat$datP[, -(c("Chr", "Start", "End", "PeakID"))]

	#Get gene start positions, accounting for reverse strand genes
	g.start <- dat$datE$start
	g.start[dat$datE$strand == "-"] <- dat$datE$end[dat$datE$strand == "-"]

	#️Create search spaces around gene start positions
	tx.s <- as.integer(g.start - txWidth/2) 
	tx.s[tx.s < 1] <- 1
	tx.e <- as.integer(g.start + txWidth/2)

	#Split gene expression dataframe into gene search space information and gene expression values dataframes
	datEa <- data.table(chr=dat$datE$chr, g.start, tx.s, tx.e, strand=dat$datE$strand, symbol=dat$datE$symbol, stringsAsFactors=F)
	datEd <- dat$datE[, -(1:5)]

	#Ordering the gene search space information and gene expression values dataframes
	tmp.order <- order(datEa$g.start)
	datEa <- datEa[tmp.order, ]
	datEd <- datEd[tmp.order, ]
	tmp.order <- str_order(datEa$chr, numeric=T)
	datEa <- datEa[tmp.order, ]
	datEd <- datEd[tmp.order, ]

	#Log2 transforming the gene expression values
	datEd <- log2(datEd + 1)

	#Calculate gene expression and peak values variance
	datE.var <- apply(datEd, 1, var)
	datP.var <- apply(datPd, 1, var)

	#Calculate gene expression and peak calling x% quantiles
	min.var.E <- quantile(datE.var, min.var.quantile.E) 
	min.var.P <- quantile(datP.var, min.var.quantile.P)

	#Remove x% lowest entries from the dataframes based on variance
	datEa2 <- datEa[which(datE.var > min.var.E), ]
	datEd2 <- datEd[which(datE.var > min.var.E), ]
	datPa2 <- datPa[which(datP.var > min.var.P), ]
	datPd2 <- datPd[which(datP.var > min.var.P), ]
        
	#Change gene search space dataframe into Genomic Ranges object
	datEa2.gr <- GRanges(seqnames=datEa2$chr, ranges=IRanges(start=datEa2$tx.s, end=datEa2$tx.e))
	names(datEa2.gr) <- datEa2$symbol
	seqlevels(datEa2.gr) <- chrs

	#Change peak information dataframe into Genomic Ranges object
	datPa2.gr <- GRanges(seqnames=datPa2$Chr, ranges=IRanges(start=datPa2$Start, end=datPa2$End))
	names(datPa2.gr) <- datPa2$PeakID
	seqlevels(datPa2.gr) <- chrs

	#Find overlaps between the gene search spaces and peaks
	fo2 <- findOverlaps(datEa2.gr, datPa2.gr)
	
	if (length(fo2) > 0) {

		#Get unique query (gene search space) and subject (peak) hits
		idxE2 <- unique(queryHits(fo2))
		idxP2 <- unique(subjectHits(fo2))

		#Subset dataframes with hits only
		datEa3 <- datEa2[idxE2, ]
		datPa3 <- datPa2[idxP2, ]
		datEd3 <- t(datEd2[idxE2, ])
		datPd3 <- t(datPd2[idxP2, ])

		#Subset Genomic Ranges object with hits only
		datEa3.gr <- datEa2.gr[idxE2]
		datPa3.gr <- datPa2.gr[idxP2]

		#Find overlaps between unique hits only
		fo3 <- findOverlaps(datEa3.gr, datPa3.gr)

		#Calculate cumulative number of overlapping gene search space hits per chromosome
		E.chr.rle <- rle(datEa3$chr)
		idx1E <- cumsum(E.chr.rle$lengths)
		idx0E <- c(1L, idx1E[-length(idx1E)] + 1L)

		#Calculate cumulative number of overlapping peak hits per chromosome
		P.chr.rle <- rle(datPa3$Chr)
		idx1P <- cumsum(P.chr.rle$lengths)
		idx0P <- c(1, idx1P[-length(idx1P)] + 1)

		return(list(datEa=datEa3, datEd=datEd3, datPa=datPa3, datPd=datPd3, idxE=queryHits(fo3), idxP=subjectHits(fo3), idx0E=idx0E, idx1E=idx1E, idx0P=idx0P, idx1P=idx1P, seed=seed))
	} else {
		return(list(datEa=NULL, datEd=NULL, datPa=NULL, datPd=NULL, idxE=NULL, idxP=NULL, idx0E=NULL, idx1E=NULL, idx0P=NULL, idx1P=NULL, seed=NULL))
	}
}


h5mat <- function(fname) {
	#' Load .h5 data as matrix
	#'
	#' @param fname character .h5 filename
	#'
	#' @return Matrix of peak-gene link correlations/samples

	#Load .h5 data
	x <- h5dump(fname)

	#Create output matrix
	res <- matrix(NA, nrow=length(x[[1]]), ncol=length(x))
	for (i in 1:length(x)) {
			res[, i] <- x[[i]]
	}
	colnames(res) <- names(x)
	res <- res[, str_order(names(x), numeric=T)]
	colnames(res) <- NULL

	return(res)
}


#execJULIA <- function(dat, N=10000L, mc.cores=NULL, work_dir=".", saveRndIdx=F, verbose=T, fname_julia="/opt/P2GL/peakToGeneLinks.jl") {

execJULIA <- function(dat, N=10000L, mc.cores=NULL, work_dir=".", saveRndIdx=F, verbose=T, fname_julia=file.path(getwd(), "peakToGeneLinks.jl")) {
	#' Calculate (spurious) peak-gene link correlations
	#'
	#' @param dat list List of filtered and ordered gene expression and peak calling dataframes
	#' @param N integer Gene start search space range
	#' @param mc.cores double Bottom percentage of genes to remove from analysis
	#' @param work_dir double Bottom percentage of peaks to remove from analysis
	#' @param saveRndIdx integer Random seed number
	#' @param verbose logical Verbose messages boolean
	#' @param fname_julia character Path to Julia script
	#'
	#' @return List of gene search space and peak links dataframes and Genomic Ranges objects

	#Initialising variables
	random.txt <- paste(unclass(as.POSIXct(Sys.time())), ".", Sys.getpid(), sep="")
	fname_tmp_RData <- paste(work_dir, "/R2J", random.txt, ".RData", sep="")
	fname_tmp_cor <- paste(work_dir, "/cor.", random.txt, ".h5", sep="")
	fname_tmp_nulldist <- paste(work_dir, "/nulldist.", random.txt, ".h5", sep="")
	fname_tmp_raw <- paste(work_dir, "/raw.", random.txt, ".h5", sep="")

	#Save list elements as seperate objects, as the save command can only save 'whole' objects
	datEa <- dat$datEa
	datEd <- dat$datEd
	datPa <- dat$datPa
	datPd <- dat$datPd
	idxE <- as.integer(dat$idxE)
	idxP <- as.integer(dat$idxP)
	idx0E <- as.integer(dat$idx0E)
	idx1E <- as.integer(dat$idx1E)
	idx0P <- as.integer(dat$idx0P)
	idx1P <- as.integer(dat$idx1P)
	seed <- as.integer(dat$seed)
	N <- as.integer(N)

	if (is.null(mc.cores)) {
		mc.cores <- detectCores()
	}

	if (saveRndIdx) {
		fname_tmp_rndIdx <- paste(work_dir, "/rndIdx.", random.txt, ".h5", sep="")
	} else {
		fname_tmp_rndIdx <- ""
	}

	if (verbose) {
		cat("transferring data (R -> Julia)", "\n")
	}

	#Saving the R dataframes and Genomic Ranges objects to RData file
    save(datEa, datEd, datPa, datPd, idxE, idxP, idx0E, idx1E, idx0P, idx1P, N, fname_tmp_cor, fname_tmp_nulldist, fname_tmp_raw, fname_tmp_rndIdx, seed, verbose, file=fname_tmp_RData)
    Sys.setenv("JULIA_NUM_THREADS"=mc.cores)

	if (verbose) {
		cat("calculating null distribution", "\n")	
		cat("mc.cores = " , mc.cores, "\n")
	}	
	
	#Using Julia to calculate (spurious) peak-gene links
	system(paste("julia ", fname_julia, " --data ", fname_tmp_RData, sep=""))

	return(list(tmp_cor=fname_tmp_cor, tmp_null=fname_tmp_nulldist, tmp_raw=fname_tmp_raw, tmp_rnd=fname_tmp_rndIdx))
}


calc_pval <- function(cor_E_P, nulldist, dat, alternative="two.sided", verbose=T) {
	#' Calculate statistical significance of peak-gene link correlations
	#' 
	#' @param cor_E_P double Vector of correlations of peak-gene links within a 0.5Mbp range
	#' @param nulldist list List of spurious peak-gene link correlation vectors
	#' @param dat list List of all peak-gene links
	#' @param alternative character Option fo which p-values to consider
	#' @param verbose logical Verbose messages boolean
	#' 
	#' @return Dataframe summarising peak-gene links correlation data 

	#Get the mean and std correlation values between gene and random peaks
	null.mean <- as.vector(nulldist$mean)
	null.sd <- as.vector(nulldist$std)

	null.mean <- null.mean[dat$idxE]
	null.sd <- null.sd[dat$idxE]


	if (verbose) {
		cat("calculating p-value & FDR", "\n")
	}

	#Calculate p-value between peak-gene links in a 0.5Mbp range and spurious peak-gene links
	if (alternative == "two.sided") {
		pval <- (1 - pnorm(abs(cor_E_P - null.mean) / null.sd, 0, 1)) * 2
	} else if(alternative == "greater") {
		pval <- (1 - pnorm((cor_E_P - null.mean) / null.sd, 0, 1))
	} else if(alternative == "less") {
		pval <- pnorm((cor_E_P - null.mean) / null.sd, 0, 1)
	} else {
		pval <- (1 - pnorm(abs(cor_E_P - null.mean) / null.sd, 0, 1)) * 2
	}

	#Summarise peak-gene links correlation data into a dataframe
	res <- cbind(dat$datEa[dat$idxE, ], idxE=dat$idxE,  dat$datPa[dat$idxP, c("Start", "End", "PeakID")], idxP=dat$idxP, r=cor_E_P, null_mean=null.mean, null_sd=null.sd, p_value=pval, FDR=p.adjust(pval, "fdr"), stringsAsFactors=F)
	res <- res[order(res$p_value), ]
	return(res)
}


findIndex <- function(dat0, dat, verbose=T) {
	#' Update information dataframes with index information
	#' 
	#' @param dat0 list List containing the original gene expression and peak calling dataframes
	#' @param dat list Selected peak - gene links dataframe
	#' @param verbose logical Verbose messages boolean
	#' 
	#' @return List of updated dataframes

	if (verbose) {
		cat("finding links between input data and results", "\n")
	}

	#Initialising variables
	datEa <- dat$datEa
	datPa <- dat$datPa

	#Add structured dataframe index to gene expression dataframe
	datE <- dat0$datE[, c("symbol", "chr", "strand", "start", "end")]
	tmpE <- 1:nrow(datEa)
	names(tmpE) <- datEa$symbol
	tmpE2 <- rep(0L, nrow(datE))
	names(tmpE2) <- datE$symbol
	tmpE2[names(tmpE)] <- tmpE
	datE <- cbind(datE, idxE=tmpE2)

	#Add structured dataframe index to peak calling dataframe
	datP <- dat0$datP[, c("PeakID", "Chr", "Start", "End")]
	tmpP <- 1:nrow(datPa)
	names(tmpP) <- datPa$PeakID
	tmpP2  <- rep(0L, nrow(datP))
	names(tmpP2) <-datP$PeakID
	tmpP2[names(tmpP)] <- tmpP
	datP <- cbind(datP, idxP=tmpP2)

	res <- list(datE=datE, datP=datP)
	return(res)
}



divIdx <- function(x, len) {
	#' Calculate sheets and entries per sheet needed
	#' 
	#' @param x integer Number of data entries
	#' @param len integer Maximum number of entries allowed per sheet
	#' 
	#' @return List of start end entries per sheet

	r1 <- ceiling(x / len)
	r2 <- x %% len

	if (r1 == 1) {
		idx0 <- 1
		idx1 <- x
	} else if (r1 > 1) {
		idx1 <- seq(len, x, len)
		idx0 <- idx1 - len + 1	
		
		if (r2 > 0) {
			idx0 <- c(idx0, idx0[length(idx0)] + len)
			idx1 <- c(idx1, x)
		}
	}
	return(list(idx0=idx0, idx1=idx1))
}



saveAsExcel <- function(result, dat_E_P, fname="results.xlsx", perSheet=1000000L, cutoff=NULL, verbose=T) {
	#' Output results to Excel sheets
	#'
	#' @param result list Peak - gene link dataframe
	#' @param dat_E_P list List containing indexed gene and peak information dataframes
	#' @param fname character Output file name
	#' @param perSheet integer Maximum number of entries in a sheet
	#' @param verbose logical Verbose messages boolean
	#' 
	#' @return Excel file with results spread over multiple sheets

	if (verbose) {
		cat("creating .xlsx book", "\n")
	}

	#Initialising variables
	datE <- dat_E_P$datE
	datP <- dat_E_P$datP
	tmp <- divIdx(nrow(result), perSheet)
	idx0 <- tmp$idx0
	idx1 <- tmp$idx1
	wb <- createWorkbook()
	modifyBaseFont(wb, fontName="MS PGothic")

	if (length(idx0) > 1) {
		result.name <- paste("result.", 1:length(idx0), sep="")
	} else {
		result.name <- "result"
	}

	#Output 'perSheet' peak - gene links per sheet
	for (i in 1:length(idx0)) {
		if (verbose) {
			cat(paste(" - writing sheet ", result.name[i], sep=""), "\n")
		}
		addWorksheet(wb, result.name[i], zoom=85)
		idx.s <- which(names(wb) == result.name[i])
		writeDataTable(wb, sheet=result.name[i], x=result[idx0[i]:idx1[i], ], withFilter=T, tableStyle="none")
		rows.idx <- 1:(idx1[i]-idx0[i]+2)
		addStyle(wb, idx.s, style=createStyle(numFmt="#,#"), rows=rows.idx, cols=which(is.element(colnames(result), "g.start")), stack=T)
		addStyle(wb, idx.s, style=createStyle(numFmt="#,#"), rows=rows.idx, cols=which(is.element(colnames(result), "tx.s")), stack=T)
		addStyle(wb, idx.s, style=createStyle(numFmt="#,#"), rows=rows.idx, cols=which(is.element(colnames(result), "tx.e")), stack=T)
		addStyle(wb, idx.s, style=createStyle(numFmt="#,#"), rows=rows.idx, cols=which(is.element(colnames(result), "Start")), stack=T)
		addStyle(wb, idx.s, style=createStyle(numFmt="#,#"), rows=rows.idx, cols=which(is.element(colnames(result), "End")), stack=T)
		addStyle(wb, idx.s, style=createStyle(numFmt="TEXT"), rows=rows.idx, cols=which(is.element(colnames(result), c("symbol"))), stack=T)
		freezePane(wb, idx.s, firstActiveRow=2, firstActiveCol=which(colnames(result) == "idxP")+ 1)
		setColWidths(wb, idx.s, cols=which(is.element(colnames(result), c("g.start", "tx.s", "tx.e", "Start", "End"))), widths="auto")
	}

	#Output gene expression information
	if (verbose) {
		cat(" - writing sheet Expression", "\n")
	}
	addWorksheet(wb, "Expression", zoom=85)
	idx.s <- which(names(wb) == "Expression")
	writeDataTable(wb, sheet="Expression", x=datE, withFilter=T, tableStyle="none")
	rows.idx <- 1:(nrow(datE)+1)
	addStyle(wb, idx.s, style=createStyle(numFmt="#,#"), rows=rows.idx, cols=which(is.element(colnames(datE), "start")), stack=T)
	addStyle(wb, idx.s, style=createStyle(numFmt="#,#"), rows=rows.idx, cols=which(is.element(colnames(datE), "end")), stack=T)
	addStyle(wb, idx.s, style=createStyle(numFmt="TEXT"), rows=rows.idx, cols=which(is.element(colnames(datE), c("symbol"))), stack=T)
	setColWidths(wb, idx.s, cols=which(is.element(colnames(datE), c("start", "end"))), widths="auto")
	freezePane(wb, idx.s, firstRow=T)
 
	#Output peak calling information
	if (verbose) {
		cat(" - writing sheet Peak", "\n")
	}
	addWorksheet(wb, "Peak", zoom=85)
	idx.s <- which(names(wb) == "Peak")
	writeDataTable(wb, sheet="Peak", x=datP, withFilter=T, tableStyle="none")
	rows.idx <- 1:(nrow(datP)+1)
	addStyle(wb, idx.s, style=createStyle(numFmt="#,#"), rows=rows.idx, cols=which(is.element(colnames(datP), "Start")), stack=T)
	addStyle(wb, idx.s, style=createStyle(numFmt="#,#"), rows=rows.idx, cols=which(is.element(colnames(datP), "End")), stack=T)
	setColWidths(wb, idx.s, cols=which(is.element(colnames(datP), c("Start", "End"))), widths="auto")
	freezePane(wb, idx.s, firstRow=T)

	if (verbose) {
		cat(" - writing", basename(fname), "\n")
	}
	saveWorkbook(wb, file=fname, T)
}


plot_cor <- function(idx, dat, nulldist, result, verbose=T) {
	
	#Subset data to only contain data from specified indices
	result <- result[idx, ]
	nulldist <- nulldist[, result$idxE, drop=F]
	datE <- dat$datEd[, result$idxE, drop=F]
	datP <- dat$datPd[, result$idxP, drop=F]

	par(mfrow=c(1, 2), mar=c(3.2, 3.2, 1.5, 0.5), mgp=c(2.0, 0.7, 0), oma=c(0,0,2,0))

	#Creating the images
	for (i in 1:length(idx)) {
		if (verbose) {
			cat(" page ", i, "\r")
		}
		x <- datE[, i]
		y <- datP[, i]

		#Plotting the correlation scatterplot
#		plot(x, y, xlab="Expression", ylab="Peak", main="correlation")	
#       edit x,y
		plot(x, y, xlab="Expression (TPM)", ylab="Peak (CPM)", main="correlation")	
		abline(lm(y ~ x), col=2)
		legend("bottomright", legend=paste("r = ", formatC(result$r[i], digits=3), sep=""), bty="n")

		#Plotting the null distribution histogram
		#edit for density for freq
		hist(nulldist[, i], breaks=seq(-1, 1, 0.02), main="null distribution", xlab="r", freq=F)
#		hist(nulldist[, i], breaks=seq(-1, 1, 0.02), main="null distribution", xlab="r", freq=TRUE)
		abline(v=result$r[i], col=2)
		legend("topleft", legend=c(paste("p = ", formatC(result$p_value[i], digits=3), sep=""), paste("FDR = ", formatC(result$FDR[i], digits=3), sep="")), bty="n")
		text = paste(result$symbol[i], " vs ", result$PeakID[i], " (", result$chr[i], ": ", result$Start[i], " - ", result$End[i], ")", sep="")
		mtext(text, side=3, cex=1.2, outer=T, line=0.5)
	}
}
#!/usr/local/bin/Rscript

suppressMessages(library(optparse))


get_path <- function() {
    return(getwd())
}

#Initialising variables
path <- get_path()

#Load necessary functions
source(paste(c(path, "peakToGeneLinksFunc.R"), collapse="/"))

option_list <- list(
	make_option(c("--fnameE", "-E"), action="store", type="character", default=NULL),
	make_option(c("--fnameP", "-P"), action="store", type="character", default=NULL),
	make_option(c("--fnameO", "-O"), action="store", type="character", default="results"),
	make_option(c("--sdE", "-e"), action="store", type="numeric", default=0.25),
	make_option(c("--sdP", "-p"), action="store", type="numeric", default=0.25),
	make_option(c("--txWidth", "-w"), action="store", type="integer", default=500000L),
	make_option(c("--alternative", "-a"), action="store", type="character", default="two.sided"),
	make_option(c("--cores", "-c"), action="store", type="integer", default=NULL),
	make_option(c("--NSim", "-n"), action="store", type="integer", default=10000L),
	make_option(c("--cutoff", "-u"), action="store", type="numeric", default=NULL),
	make_option(c("--seed", "-s"), action="store", type="integer", default=NULL),
	make_option(c("--saveRndIdx", "-r"), action="store", type="logical", default=F),
	make_option(c("--verbose", "-v"), action="store", type="logical", default=T),
	make_option(c("--saveExcel", "-x"), action="store", type="logical", default=FALSE)  # デフォルトを FALSE に設定
)

#Parse command line arguments
opt <- parse_args(OptionParser(option_list=option_list))

#Create the output directory if it does not yet exist
dir.create(file.path(path, "output"), showWarnings = FALSE)
 
#Load gene expression and peak calling files as dataframes
dat0 <- load_data(opt$fnameE, opt$fnameP, opt$verbose)

#Set random seed, if not provided
if (is.null(opt$seed)) {
	seed <- floor(runif(1) * 2^30)
} else {
	seed <- opt$seed
}

print(seed)

#Find links between genes and peaks
dat <- find_link(dat=dat0, txWidth=opt$txWidth, min.var.quantile.E=opt$sdE, min.var.quantile.P=opt$sdP, seed=seed, verbose=opt$verbose)

if (is.null(dat$datEa)) {
	print("no candidate links between genes and peaks!")
	quit("no")
}

#Calculate (spurious) peak-gene link correlations
res_julia <- execJULIA(dat=dat, N=opt$NSim, mc.cores=opt$cores, work_dir=path, saveRndIdx=opt$saveRndIdx, verbose=opt$verbose, fname_julia=paste(c(path, "peakToGeneLinks.jl"), collapse="/"))

if (opt$verbose) {
	cat("transferring data (Julia -> R)", "\n")
}

#Load peak-gene correlation value files as vectors
cor_E_P <- as.vector(h5dump(res_julia$tmp_cor)$x)
nulldist <- h5dump(res_julia$tmp_null)

#Save correlation value vectors into a list
if (opt$saveRndIdx) {
	rndIdx <- h5mat(res_julia$tmp_rnd)	
}

#Calculate peak-gene link correlation p-values
result <- calc_pval(cor_E_P=cor_E_P, nulldist=nulldist, dat=dat, alternative=opt$alternative, verbose=opt$verbose)

#Remove temporary files
unlink(list.files(path, pattern="cor.|nulldist."))
unlink(list.files(path, pattern=".RData"))

#Subsetting dataframe to only contain values lower than cutoff value
if (!is.null(opt$cutoff)) {
	result <- result[result$FDR <= opt$cutoff]
	opt$fnameO <- paste(opt$fnameO, as.character(opt$cutoff), sep="_")

	if (length(result$FDR) == 0) {
		cat("no links were found for significance level", opt$cutoff, "\n")
		quit()
	}
}

#Add index to information dataframes
dat_E_P <- findIndex(dat0=dat0, dat=dat, verbose=opt$verbose)

# Save results to Excel if the option is set
if (opt$saveExcel) {
	saveAsExcel(result=result, dat_E_P=dat_E_P, fname=paste(path, "/output/", opt$fnameO, ".xlsx", sep=""), verbose=opt$verbose)
}
if (opt$verbose) {
	cat("writing", paste(opt$fnameO, ".tsv", sep=""), "\n")
}
write.table(result, file=paste(path, "/output/", opt$fnameO, ".tsv", sep=""), sep="\t", row.names=F, col.names=T, quote=F)

if (opt$verbose) {
	cat("writing", paste(opt$fnameO, ".RData", sep=""), "\n")
}
save(dat, result, dat_E_P, res_julia, file=paste(path, "/output/", opt$fnameO, ".RData", sep=""))

if (opt$verbose) {
	cat("finished!", "\n")
}

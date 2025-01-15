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
make_option(c("--fnameI", "-I"), action="store", type="character", default="results"),
make_option(c("--fnameO", "-O"), action="store", type="character", default="results"),
make_option(c("--idx", "-x"), action="store", type="character", default="1:100"),
make_option(c("--idxFile", "-f"), action="store", type="character", default=NULL),
make_option(c("--verbose", "-v"), action="store", type="logical", default=T)
)

#Parse command line arguments
opt <- parse_args(OptionParser(option_list=option_list))

#Checking input and output file names
if (!grepl(".RData$", opt$fnameI)) {
	opt$fnameI <- paste(opt$fnameI, ".RData", sep="")
}
if (!grepl(".pdf$", opt$fnameO)) {
	opt$fnameO <- paste(opt$fnameO, ".pdf", sep="")
}
if (opt$verbose) {
	cat("loading", opt$fnameI, "\n")
}

#Reload saved data
load(opt$fnameI)

#Create index vector or load file of links to visualise
if (is.null(opt$idxFile)) {
	idx <- eval(parse(text=opt$idx))
} else {
	idx <- scan(opt$idxFile, integer(0))
}

if (opt$verbose) {
	cat("writing", opt$fnameO, "\n")
}

#Create pdf containing images
pdf(file=paste(path, "/output/", opt$fnameO, sep=""), width=11.7, height=8.3, onefile=T)
plot_cor(idx, dat, h5mat(res_julia$tmp_raw), result, opt$verbose)
dev.off()

#Remove temporary files
unlink(list.files(path, pattern="raw."))

if (opt$verbose) {
	cat("finished!", "\n")
}


#!/usr/local/bin/julia

using ArgParse
using Base.Threads
using HDF5
using Random
using RData
using StatsBase
using Statistics


"""
Calculate correlation between peak - gene links

### Input
 - `datE`  -- Matrix of gene expression values
 - `datP`  -- Matrix of peak values
 - `fname` -- Peak-gene links correlation output filename

### Output
 .h5 file containing the correlation values between peak-gene links
"""
function calcCor(datE, datP, fname="")

	#Initialising variables
	numE = size(datE)[2]
	res = Array{Float64}(undef, numE)

	#Calculate correlation between all peak-gene links
	@threads for i = 1:numE
		res[i] = cor(view(datE, :, i), view(datP, :, i))
	end

	#Output corelation results
	if fname != ""
		f = h5open(fname, "w");
		f["x"] = res
		close(f)
	end	
end


"""
Set random seeds for when random sampling spurious peak-gene link assocations

### Input
 - `len`  -- Number of seeds to set
 - `seed` -- Seed number
			 (default: Int.(floor.(rand(1)*2^30))[1])

### Output
Seed numbers vector
"""
function setSeeds(len, seed=Int.(floor.(rand(1)*2^30))[1])
	Random.seed!(seed)
	res = Int.(floor.(rand(len) * 2^30))
	return res
end


"""
Output spurious peak-gene link correlations/samples per chromosome

### Input
 - `x` 	   -- Peak-gene link correlations/samples results array
 - `fname` -- Output filename

### Output
Peak-gene link correlations/samples .h5 file
"""
function saveEachCol(x, fname)

	#Initialising variables
	N = size(x)[2]
	dataName = map(string, repeat("x", N), 1:N);

	#Outputting
	f = h5open(fname, "w");
	for i in 1:N
		f[dataName[i]] = x[:, i];
	end
	close(f);
end


"""
Calculate correlation between random spurious peak-gene links

### Input
 - `datE`	 -- Matrix of gene expression values
 - `datP`	 -- Matrix of peak values
 - `idx0E`	 -- Gene links per chromosome starting index vector
 - `idx1E`	 -- Gene links per chromosome ending index vector
 - `idx0P`	 -- Peak links per chromosome starting index vector
 - `idx1P`	 -- Peak links per chromosome ending index vector
 - `N`		 -- Number of random links to take
				(default: 10000)
 - `fname`   -- Correlation mean & std output filename
				(default: "")
 - `fname2`	 -- Randomly sampled peak-gene links
				(default: "")
 - `f_raw` 	 -- Correlation raw output output filename
				(default: "")
 - `seed` 	 -- Seed number
				(default: Int.(floor.(rand(1)*2^30))[1])
 - `verbose` -- Verbose messages boolean
				(default: true)

### Output
Spurious peak-gene link correlations
"""
function randomCor(datE, datP, idx0E, idx1E, idx0P, idx1P, N=10000, fname="", fname2="", f_raw="", seed=setSeeds(1)[1], verbose=true)

	#Initialising variables
	seeds = setSeeds(size(datE)[2], seed)
	numE = size(datE)[2]
	res = Array{Float64}(undef, N, numE)

	if fname2 != ""
		rndIdx = Array{Int32}(undef, N, numE)
	end

	if verbose
		println("chromosome")
	end

	#Loop over every chromosome
	for i = 1:length(idx0E)
		if verbose
			println(" - " * string(i))
		end

		#Select all genes on the chromosome and peaks not on the chromosome
		idxE = idx0E[i]:idx1E[i]
		idxP = setdiff(1:size(datP)[2], idx0P[i]:idx1P[i])

		#Calculate correlation between N randomly sampled spurious peak-gene link associations
		@threads for j in idxE
			Random.seed!(seeds[j])
			idxP_rnd = sample(idxP, N, replace=false)
			res0 = Array{Float64}(undef, N)
			for k = 1:N
				res0[k] = cor(view(datE, :, j), view(datP, :, idxP_rnd[k])) 
			end
			res[:, j] = res0
			if fname2 != ""
				rndIdx[:, j] = idxP_rnd 
			end
		end
	end

	#Output randomly sampled spurious peak-gene links
	if fname2 != ""
		saveEachCol(rndIdx, fname2)
	end

	#Output spurious correlation's mean & std results
	if fname != ""
		f = h5open(fname, "w");
		f["mean"] = mean(res, dims=1);
		f["std"] = std(res, dims=1);
		close(f);
	end

	#Output raw spurious correlation results
	if f_raw != ""
		saveEachCol(res, f_raw)
	end
end


"""
Parse command line arguments

### Output
Returns parsed arguments
"""
function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin
		"-d", "--data"
			help = "Peak - gene search space links Rdata file"
			arg_type = String
			required = true
	end
	return parse_args(s)
end


#Parse command line arguments
parsed_args = parse_commandline()

#Loading R data
objs = load(parsed_args["data"]);

#Calculate (spurious) peak-gene links
calcCor(objs["datEd"][:, objs["idxE"]], objs["datPd"][:, objs["idxP"]], objs["fname_tmp_cor"])
randomCor(objs["datEd"], objs["datPd"], objs["idx0E"], objs["idx1E"], objs["idx0P"], objs["idx1P"], objs["N"], objs["fname_tmp_nulldist"], objs["fname_tmp_rndIdx"], objs["fname_tmp_raw"], objs["seed"], objs["verbose"])

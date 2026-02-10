#setwd("../")
library(ape)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(gplots)
library(ggridges)
library(ggpubr)
library(coda)
library(lubridate)
library(scales)
library(latex2exp)
library(pbapply)
#library(bdskytools) # Installed from GitHub. Requires devtools
library(svglite)
library(cowplot)
library(grid)
library(emmeans)

# For the phylogeo
library(igraph)
library(sf)
library(sfnetworks)
library(tidygraph)
library(GGally)
#library(network)
library(netmap)
library(ggraph)
#library(sna)

# For Ali stuff only
#library(fitdistrplus)

setwd("../")
# Functions
kneedle <- function(x,guess=length(x)){
       	# This is an implementation by Zachery Dickson. Comes from https://raghavan.usc.edu/papers/kneedle-simplex11.pdf 
    if(is.na(guess) || guess > length(x) / 2){
        guess = length(x)/2
    }
    x = x[1:(guess*2)]
    y = lowess(x)$y
    y0 = y[1]
    yn= y[length(y)]
    x0=0
    xn=length(y)-1
    m = (yn-y0)/(xn-x0)
    b = y0
    d = abs(y - (m*(x0:xn) + b)) * asin(pi/4)
    return(which.max(d) + 1)
}

beastLogPartitioning <- function(logFile, partitionNames){ # For when trying to deal with all the possible partitions

	if(length(partitionNames) == 1){ # If there's only one partition!
		chainList <- list(logFile)
	}else{
		partitions <- lapply(colnames(logFile), function(x){# This is going to help me separate the two partitions. IMPORTANT the first three are columns are going to be the same for both!
			y <- strsplit(x, split = "\\.") |> unlist()
			return(y[2])
					})
		sharedColumns <- which(is.na(partitions))
		
		# Now we're getting the columns split into the requisite number of groups
		columnPartition <- unlist(partitions) |> factor()
		columnPartition <- split(seq_along(columnPartition), columnPartition) # The + 3 is because three of the variables are shared
		
		chainList <- list()
		for(i in 1:length(columnPartition)){ # This is going to be the final part of getting it to work
			tmp <- logFile[,c(sharedColumns, columnPartition[[i]])] # Splitting the DF
			colnames(tmp) <- colnames(tmp) |> sub(pattern = "\\.\\d",replacement = "")
			chainList[[i]] <- tmp
		}
		rm(tmp,i)

	}

	names(chainList) <- partitionNames	
	
	return(chainList)
}

BeastLogAnalysis <- function(files, partitionNames, burnInPeriod = NA, diagPlotting = T, plotting = T){
	if(!any(c(exists("files"),exists("protein")))){
		stop("Missing an argument!")
	}
	dir.create(path = "Results", showWarnings = F) # Only care if the folder doesn't exist. Don't warn otherwise
	dir.create(path = "Figures", showWarnings = F)

	clusterCount <- min(length(files), parallel:::detectCores() - 1) # How many cores should I use?

	cat("Reading the Beast Files\n")
	logFiles <- pblapply(files, cl = clusterCount,function(x){
				     df <- read.delim(x, comment.char = "#")
				     if(colnames(df)[1] == "state"){
					colnames(df)[1] <- "Sample"
				     }
					colnames(df) <- gsub("\\.$","",colnames(df))

				     return(df)
					})
	
	# Calculating some smaller things
	maxSize <- logFiles[[1]]$Sample |> max()
	stateSize <- logFiles[[1]]$Sample[3] - logFiles[[1]]$Sample[2]

	# Figuring out the burn-in
	if(is.na(burnInPeriod)){
		cat("Figuring out the burn-in with Kneedle\n")
		guessSize <- logFiles[[1]]$Sample[nrow(logFiles[[1]])] / stateSize
		
		burnInPeriod <- max(pbsapply(logFiles, cl = clusterCount, function(x){kneedle(x$joint, guess = guessSize) * stateSize})) # We want to choose the largest here
	# to make sure it's applicable to all our chains
	}else{
		if(burnInPeriod >= max(logFiles[[1]]$Sample)){
		       	errorCondition("The burnInPeriod is larger than the log size!")
		}
	}
	
	# Let's filter it down so that we remove the burn-in!
	beastLogs <- pblapply(logFiles, cl = clusterCount, function(x){
				      burninIndex <- max(which(x$Sample < burnInPeriod)) + 1
				      return(x[burninIndex:nrow(x),])
				}) # the -1 here is removing the Samples column

	# Need to remove the monotonic sections as they're not interesting. Should be consistent across all runs so only looking at first
	ind <- sapply(beastLogs[[1]], function(x){
				    y <- ifelse(mean(x) == 0 | is.infinite(mean(x)) | is.na(mean(x)) | any(abs(x) > 1e100 | abs(x) < 1e-100), F, T)
				    return(y)
				})
	beastLogs <- lapply(beastLogs, function(x){return(x[ind])})
	
	cat("Converting the tibbles to mcmc. Also concatenating the chains\n")
	beastMCMCs <- pblapply(beastLogs, cl = clusterCount, function(x) x[,-1] |>
		  	        mcmc(start = burnInPeriod, end = max(logFiles[[1]]$Sample),thin = stateSize)) #|>
				#as.mcmc.list()
	combinedChain <- beastLogs |> bind_rows() |> select(-Sample) |> mcmc(thin = stateSize)
	combinedLogs <- combinedChain |> as_tibble()

	# It is at this point where I want to start separating the chains
	cat("Separating the partitions into a list of dataframes\n")
	logList <- beastLogPartitioning(combinedLogs, partitionNames)
	chainList <- beastLogPartitioning(combinedChain, partitionNames)

	#combinedChain <- append_chains(beastMCMCs |> unlist())
	cat("Calculating the ESS and other standard summary statistics\n")
	beastSummarizedList <- pblapply(chainList, cl = clusterCount, function(part){
		tmp <- tibble(Variable = colnames(part),
				  ESS = effectiveSize(part),
				  Mean = part |> as_tibble() |> sapply(mean),
				  SD = part |> as_tibble() |> sapply(sd), 
				  Median = part |> as_tibble() |> sapply(median), 
				  HPDLo = HPDinterval(part)[,1], 
				  HPDHi = HPDinterval(part)[,2])
		return(tmp)
		})

	 # This is going to print each partition separately. Easier this way to know what you're looking at. Going to combine this with the low ESS filter
	 # since I'd have to make an apply command there too
	for(protein in partitionNames){
		df <- beastSummarizedList[[protein]]
		write.table(df, paste0("Results/",protein,"_Summary.tab"), sep = "\t", quote = F, row.names =F) # Printing everything

		if(any(df$ESS < 200)){
			lowESS <- df$Variable[which(df$ESS < 200)]
			warn <- paste0("The following variables had an ESS less than 200: ", paste(lowESS, collapse = ", "), ". Please check the summary and trace files.")
			warning(warn)
		}
	}	
	rm(df)

	if(diagPlotting){
		cat("Plotting the traces of your chains. Can take a while...\n")
		multiSummarized <- pblapply(1:length(beastLogs), cl = clusterCount, function(x){
				y <- beastMCMCs[[x]]
				out <- 	tibble(Variable = colnames(y),
						  ESS = effectiveSize(y),
						  Mean = y |> as_tibble() |> sapply(mean),
						  SD = y |> as_tibble() |> sapply(sd), 
						  Median = y |> as_tibble() |> sapply(median), 
						  HPDLo = HPDinterval(y)[,1], 
						  HPDHi = HPDinterval(y)[,2],
						  Chain = x)
				return(out)
					  }) |> bind_rows() #|> t() |> as_tibble()
		# Now we need to separate variables again. This time its going to be rows and the chains need to be pulled apart and recombined again
		multiList <- split.data.frame(multiSummarized, f = multiSummarized$Chain)

		multiList <- lapply(multiList, function(x){
			y <- x[,-1] |> t() |> as.data.frame()
			colnames(y) <- x[,1] |> pull()

			z <- beastLogPartitioning(y, partitionNames) 
			return(z)
			 })

		multiSummarizedPartitioned <- list()
		for(i in partitionNames){
			tmp = NULL
			for(j in 1:length(multiList)){
				tmp <- tmp |> bind_rows(multiList[[j]][[i]] |> t() |> as.data.frame())
			}
			tmp <- tmp |> mutate(Variable = rownames(tmp), Variable = gsub("\\.\\..*", "", Variable))
			multiSummarizedPartitioned[[i]] <- tmp
		}
		rm(tmp)

		# This is making the caterpillar plots
		diagData <- lapply(1:length(beastLogs), function(x){beastLogs[[x]] |> mutate(Chain = x)}) |> bind_rows()
		diagDataList <- beastLogPartitioning(diagData, partitionNames)
		
		pages <- ceiling(ncol(combinedLogs)/20)
		# Need to make sure I get the right size for these plots
		for(protein in partitionNames){
			cat(paste0("Trace Plot: ", protein, "\n"))
			tracePlot <- diagDataList[[protein]] |> pivot_longer(-c(Sample, Chain), names_to = "Variable") |> 
				ggplot(aes(x = Sample, y = value, colour = as.factor(Chain))) +
				scale_colour_manual(values = coloursChain, "Chain") +
				theme_classic() +
				geom_line(alpha = 0.5) +
				geom_smooth(se = F) +
				theme(legend.position = "bottom")

			pdf(paste0("Figures/",protein,"_Trace.pdf"), width = 12, height = 9)
				for(j in 1:pages){			
					print(tracePlot + facet_wrap_paginate("Variable", scales = "free_y", nrow = 4, ncol = 5, page = j))
				}
			dev.off()

#			ggsave(paste0("Figures/",protein,"_Trace.pdf"), tracePlot, width = 12, height = 9)
	
		# This is making the distribution plots
			cat(paste0("Distribution Plot: ", protein, "\n"))
			distPlot <- diagDataList[[protein]] |> pivot_longer(-c(Sample, Chain), names_to = "Variable") |>
	       			ggplot(aes(x = value, y = Chain, fill = as.factor(Chain))) +
				#facet_wrap("Variable", scales = "free") +
				scale_fill_manual(values = coloursChain, "Chain") +
				scale_colour_manual(values = coloursChain, "Chain") +
				theme_classic() +
				geom_segment(multiSummarizedPartitioned[[protein]], inherit.aes = F,
					  mapping = aes(colour = as.factor(Chain),x= HPDLo, xend = HPDHi , y = Chain)) +
				geom_point(multiSummarizedPartitioned[[protein]], inherit.aes = F,
					  mapping = aes(colour = as.factor(Chain), x= Median , y = Chain)) +
				geom_density_ridges(alpha = 0.5, colour = NA) +
				theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

			pdf(paste0("Figures/",protein,"_Density.pdf"), width = 12, height = 9)
				for(j in 1:pages){			
					print(distPlot + facet_wrap_paginate("Variable", scales = "free", nrow = 4, ncol = 5, page = j))
				}
			dev.off()
			#ggsave(paste0("Figures/",protein,"Density.pdf"), distPlot, width = 12, height = 9)
		}
	}

	if(plotting){
		cat("Preparing the Phylogeny\n")
		cat("Currently limited to using log combiner and tree annotator\n")
		# Getting the files and folders ready
		#It's a little more complicated here since we now need to find *all* the trees present for the partitions
		fileStem <- basename(gsub(".log.gz", "", files)) |> unique()
		fileTrees <- list.files(path, pattern = paste0(".*",fileStem,".*trees"), recursive = T, full.names = T)
		
		# Now to split them into their partitions
#		fileTrees <- lapply(partitionNames, function(x){
#			       return(fileTrees[grepl(x, fileTrees)])
#		})
#		names(fileTrees) <- partitionNames
		

		tmp <- strsplit(fileTrees[[1]][[1]], "/")[[1]]
		combinedTreesFolder <- paste(paste(tmp[1:(length(tmp) - 2)], collapse = "/"), "Trees/", sep = "/")
		dir.create(path = combinedTreesFolder, showWarnings = F)

		# Getting the static variables ready for the logcombiner command
		startLogCombiner <- "logcombiner -log"
		resampleSize <- stateSize * 10
		logCombinerBurnIn <- ceiling(burnInPeriod/maxSize * 100) # Due to how they wrote the program, we need to round *up* to the nearest whole integer

		# Similar story for the treeannotator
		startTreeAnnotator <- "treeannotator -burnin 0 -file" # Note that the burnin is zero since the logcombiner will take care of that for us in our trees

		cat(paste0("Creating the ", partitionNames, " tree\n"))
		# Getting the log output ready
		fileOutLogCombiner <- paste0(combinedTreesFolder, partitionNames, "Combined.trees") 
		
		# Now we need to get the files ready
		partFiles <- fileTrees
		logCombinerCommand <- paste(startLogCombiner, paste(partFiles, collapse = " -log "), "-b", logCombinerBurnIn, "-o", fileOutLogCombiner, "-resample", resampleSize)

		# Running the command
		system(logCombinerCommand)

		# Extracting the MCC
		fileOutAnnotator <- paste0(combinedTreesFolder, partitionNames, ".nexus") 
		treeAnnotatorCommand <- paste(startTreeAnnotator, fileOutLogCombiner, fileOutAnnotator)

		# Running the command
		system(treeAnnotatorCommand)
		
		# Deleting the large trees file
		rmCommand <- paste("rm -f", fileOutLogCombiner)
		system(rmCommand)
		
	}

	cat(paste0(rep("#",40), collapse = ""), "\n")
#	cat("Number of Intervals = ", intervals,"\n")
	cat("Burn-In = ", burnInPeriod,"\n")
	cat(paste0(rep("#",40), collapse = ""), "\n")
	cat("Diagnostic Plots = ", diagPlotting,"\n")
	cat("Extracted Phylogeny = ", plotting,"\n")
	cat(paste0(rep("#",40), collapse = ""), "\n")
	cat(paste0(rep("#",40), collapse = ""), "\n")
	return(list("LogFiles" = logList, "MCMCChain" = chainList, "Summarized" = beastSummarizedList))
}

bdskyReCalculator <- function(beastOutSummarized, mrsd){
	if(!any(c(exists("beastOut"),exists("mrsd")))){
		stop("Missing an argument!")
	}

	# Let's pull out the required rows from the beast summary table
	RequiredReValues <- beastOutSummarized |> filter(grepl("Tree.height|reproductiveNumber", Variable)) 

	# What's the median height of the tree?
	height <- RequiredReValues$Median[1]

	# We need to get the number of intervals for the height to be broken up into
	numIntervals <- nrow(RequiredReValues) - 1 # -1 since the tree height is included
	intervals <- height/numIntervals 

	# Calculating the dates where each dimension starts
	Dates <- rev(date_decimal(mrsd - intervals * 1:numIntervals)) # reversed so that the intervals are in the correct order
	
	# Making the dataframe!
	BDSKYReResults <- RequiredReValues[-1,] |> mutate(Start = Dates, End = c(Dates[-1],date_decimal(mrsd)))

	return(BDSKYReResults)
}

bdmmReParsing <- function(bdmmPrimeResults, mrsd){

	# Finding the median of the origin
	originbdmmPrime <- bdmmPrimeResults |> filter(grepl("origin", Variable)) |> pull(Median)

	# We need to get the dates sorted out first
	Dates <- bdmmPrimeResults |> filter(grepl("ReSP.*_endtime", Variable)) |> select(Variable, Median) |> rename(Date = Median) |>
		mutate(Date = mrsd - (originbdmmPrime - Date)) |> mutate(End = as.Date(date_decimal(Date)))
	Dates$Start <- c(mrsd - originbdmmPrime, Dates$Date[1:(nrow(Dates)-1)])
	Dates <- Dates |> mutate(Start = as.Date(date_decimal(Start))) |> mutate(Variable = gsub("_endtime","", Variable))

	# This here will make sure we have the final era included
	lastEra <- gsub(".*i","",Dates$Variable) |> as.numeric() |> max() +1
	Dates[nrow(Dates) + 1,] <- list(Variable = paste0("ReSPEpi.i", lastEra), Date = mrsd,
					End = as.Date(date_decimal(mrsd)), Start = Dates[nrow(Dates),"End"] |> pull())
	
	# Now we need to merger the Re values with the dates that we just pulled out
	reResults <- bdmmPrimeResults |> filter(grepl("ReSPEpi.i\\d*$", Variable)) |> select(Variable, Median, HPDLo, HPDHi)
	reResults <- inner_join(Dates, reResults)

	return(reResults)
}

dateConverterOntario <- function(x){
	dateSplit <- strsplit(x, " ") |> unlist()
	dateSplit[1] <- which(month.name == dateSplit[1])
	return(paste(dateSplit[c(3,1,2)], sep = "-", collapse = "-"))
}

colourProvince <- c("Canada" = "#F59499", "BC" = "#FEE1D9", "SK" = "#F1E3F2", "MB" = "#DBE9F5", "ON" = "#D1EFD4", "NB" = "#F8E5C3", "NS" = "#CFEDED")
colourProvinceLine <- c("Canada" = "#EB2D37", "BC" = "#C64A1C", "SK" = "#92278F", "MB" = "#0369AC", "ON" = "#2B8737", "NB" = "#8A600D", "NS" = "#367A76")
coloursChain <- c("Black", "#3193CC", "#F0454B","#39B54A", "#B975B7", "#F15A22", "#CBA52E")

provinceConversionList <- list("British_Columbia" = "BC", "New_Brunswick" = "NB", "Ontario" = "ON", "Manitoba" = "MB", "Nova_Scotia" = "NS", "Saskatchewan" = "SK")

#### Metadata ####
alignmentFiles <- list.files("Beast/CanadianKCladeTempestFiltered/", pattern = "*fasta", recursive = T, full.names = T)
dates <- read.delim("Beast/All_meta_v2.tsv", header = T, col.names = c("Genome", "Date", "Location")) |> as_tibble()
nextCladeData <- read.delim("Nextclade/GisaidNextclade.tab", header = T) |> select(seqName, clade) |> filter(clade == "K")

seqNames <- lapply(alignmentFiles, function(f){
	seqNames <- read.FASTA(f) |> names()
	fileName <- gsub(".*/|.fasta","",f)
	data.frame(Genome = seqNames)#, Partition = fileName)
	}) |> bind_rows()
datesSeqs <- seqNames |> left_join(dates) |> unique()

datesMRSDs <- datesSeqs |> mutate(Date = as.Date(Date)) |> group_by(Location) |> filter(Date == max(Date)) |> select(-Genome) |>
	distinct() |> pull(Date, name = Location)

numGenomes <- datesSeqs |> group_by(Location) |> count(name = "Genomes") |> mutate(Location = provinceConversionList[Location] |> unlist())
# Preparing for a plot
allweekdays <- wday(as.Date(as.Date("2025-01-01"):as.Date("2025-12-31")))
names(allweekdays) <- as.Date(as.Date("2025-01-01"):as.Date("2025-12-31"))
startWeeks <- names(allweekdays[allweekdays == 1])

ontarioSampleCollection <- datesSeqs |> tibble() |> mutate(Date = as.Date(Date), Week = startWeeks[epiweek(Date) - 1]) |> filter(Location == "Ontario") |> group_by(Week) |>
	count(name = "Ontario") |> mutate(Week = as.Date(Week)) 

startWeeks <- startWeeks |> as.Date() 


samplingData <- datesSeqs |> tibble() |> mutate(Date = as.Date(Date), Week = startWeeks[epiweek(Date) - 1])  |> group_by(Week) |>
	count(name = "Canada") |> mutate(Week = as.Date(Week)) |> left_join(ontarioSampleCollection) |>
	mutate(Ontario = replace(Ontario, is.na(Ontario), 0), Canada = Canada - Ontario) |> rename(ON = Ontario) |>
	pivot_longer(-Week, names_to = "Province", values_to = "Genomes") |>
	ggplot(aes(x = Week, y = Genomes, fill = Province)) +
	theme_classic() +
	scale_fill_manual(values = colourProvinceLine, "") +
	scale_y_continuous(breaks = extended_breaks(n = 10)) +
	scale_x_datetime(breaks = startWeeks[startWeeks > as.Date("2025-08-07")] |> as.POSIXct(), date_labels = "%Y-%m-%d") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom") +
	geom_col() +
	xlab("")


ggsave("Figures/SamplingDistribution.png", samplingData, width = 6, height = 4)
######## BDSKY #####
# Parsing the Ontario Results
path <-"Beast/Results/BDSKY" 
#files <- list.files(path, pattern = paste0(".*HA.*log.gz"), recursive = T, full.names = T)
files <- list.files(path, pattern = paste0("*Ontario_BDSKY8Dims.log.gz"), recursive = T, full.names = T)

Ontario <- BeastLogAnalysis(files, "H3BDSKY_Ontario", diagPlotting =T, plotting = F, burnInPeriod = 6e6)

OntarioBDSKYRe <- bdskyReCalculator(Ontario$Summarized$H3BDSKY_Ontario, datesMRSDs["Ontario"] |> decimal_date()) |> mutate(Province = "ON")

# Parsing the Canadian Results
files <- list.files(path, pattern = paste0("*Canada_BDSKY_12Dim.log.gz"), recursive = T, full.names = T)

Canada <- BeastLogAnalysis(files, "H3BDSKY_Canada", diagPlotting =T, plotting = F, burnInPeriod = 6e6)

CanadaBDSKYRe <- bdskyReCalculator(Canada$Summarized$H3BDSKY_Canada, max(datesMRSDs) |> decimal_date()) |> mutate(Province = "Canada")

# Let's pull out the tMRCAs for each now
Ontario$Summarized$H3BDSKY_Ontario |> filter(Variable == "Tree.height") |> mutate(across(c(Mean, Median, HPDLo, HPDHi), ~ (datesMRSDs["Ontario"] |> decimal_date()) - . )) |>
	mutate(across(c(Mean, Median, HPDLo, HPDHi), ~ date_decimal(.)))

Canada$Summarized$H3BDSKY_Canada |> filter(Variable == "Tree.height") |> mutate(across(c(Mean, Median, HPDLo, HPDHi), ~ (max(datesMRSDs) |> decimal_date()) - . )) |>
	mutate(across(c(Mean, Median, HPDLo, HPDHi), ~ date_decimal(.)))
# Plotting the results
finalSegment <- OntarioBDSKYRe |> bind_rows(CanadaBDSKYRe) |> group_by(Province) |> filter(Start == max(Start))

OntarioBDSKYRe |> bind_rows(CanadaBDSKYRe) |> write.table(file = "Results/BDSKYResults.tab", sep = "\t", row.names = F, quote = F)

H3N2RePlot <- OntarioBDSKYRe |> bind_rows(CanadaBDSKYRe) |>
	ggplot(aes(x = Start, xend = End, xmin = Start, xmax = End, y = Median, ymin = HPDLo, ymax = HPDHi, colour=  Province, fill = Province)) +
	theme_classic() +
	scale_x_datetime(breaks = breaks_pretty(n = 6), date_labels = "%b %Y") +
	scale_y_continuous(breaks = extended_breaks(n = 5)) +
	scale_colour_manual(values = colourProvinceLine, "Area") +
	scale_fill_manual(values = colourProvince, "Area") +
	geom_rect(inherit.aes = F, aes(xmin = as.Date("2025-10-14"), xmax = as.Date("2025-10-27"), ymin = -Inf, ymax = Inf), fill = "grey70", colour = NA) +
	#geom_vline(xintercept = as.Date("2025-10-27"), lty =2, colour = "grey10") +
	geom_hline(yintercept = 1, lty =2, colour = "grey10") +
	geom_rect(alpha = 0.5, lwd = 0.25, lty = 3) +
	geom_segment(data = finalSegment) +
	geom_step() +
	annotate(geom = "text", x = as_datetime("2025-10-20 12:00:00"), y = 3.75, label = "Vaccines Released") +
	ylab(TeX("$R_t$")) +
	xlab("")  +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom", legend.title.position = "top")
ggarrange(H3N2RePlot, samplingData, ncol = 1, align = "hv", labels = "AUTO", common.legend = T, legend = "bottom")
ggsave("Figures/BDSKY_HA.png",  width = 6, height = 8)

##### Phylogeography ####
# Reading the shapefile
#canadaMap <- read_sf("ShapeFiles/lpr_000b21a_e.shp") <- # NOTE: Shapefiles not included in repository. Please see https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?year=21
meanLocations <- read.delim("ShapeFiles/Locations.tab", header = F, col.names = c("Province", "lat", "lon"))

# Reading the asymmetric Results
asymmetricResults <- read.delim("PhyloGeographicResults/AsymmetricRates.json.txt") |> 
	tibble() |> mutate(LogBayes = log10(BAYES_FACTOR)) |>
	mutate(FROM = provinceConversionList[FROM] |> unlist(), TO = provinceConversionList[TO] |> unlist()) |>
	filter(LogBayes >= 1)

provincesPresent <- unique(c(asymmetricResults$FROM, asymmetricResults$TO))
meanLocations <- meanLocations[meanLocations$Province %in% provincesPresent,]
meanLocations <- meanLocations |> left_join(numGenomes, by = c("Province" = "Location"))

# Making the network
net <- graph_from_data_frame(asymmetricResults[,c(1,2,5)])

# Reorganizing meanLocations dataset
meanLocations <- meanLocations[c(1,3,5,4,2,6),]
#net <- set_vertex_attr(net, name = "lat", index = V(net), value = meanLocations$lat)
#net <- set_vertex_attr(net, name = "lon", index = V(net), value = meanLocations$lon)
net <- set_vertex_attr(net, name = "Genomes", index = V(net), value = meanLocations$Genomes)

#
canadaMap$PREABBR <- c("NL", "PEI", "NS", "NB", "QC", "ON", "MB", "SK", "AB", "BC", "YT", "NWT", "NU")

## Trying to make the plots
legendColours <- legendGrob(c("Canada", "British Columbia", "Saskatchewan", "Manitoba", "Ontario", "New Brunswick", "Nova Scotia"), hgap = 0.75, vgap = 0.75, pch = 22, 
	nrow = 3, ncol = 3, gp = gpar(col = colourProvinceLine, fill = colourProvince, cex = 0.75))

asymNetwork <- ggraph(net, layout = "manual", x = meanLocations$lon, y = meanLocations$lat) +
	theme_bw() +
	geom_sf(inherit.aes = F, data = canadaMap, aes(fill = PREABBR)) +
	geom_edge_link2(aes(width = LogBayes),start_cap = circle(0.25,"cm"), end_cap = circle(0.25,"cm"),
			arrow = arrow(type = "closed", length = unit(0.25, "cm"))) +
	geom_node_point(aes(size = Genomes, colour = name), position = "jitter") +
	scale_edge_width(range = c(0.5,2), TeX("$\\log_{10}(BF)$")) +
	scale_size(range = c(1, 3)) +
	scale_colour_manual(values = colourProvinceLine, na.value = "white", "Province") +
	scale_fill_manual(values = colourProvince, na.value = "white", "Province") +
	guides(fill = guide_none(), colour = guide_none(), custom = guide_custom(legendColours, title = "Jurisdiction"), edge_width = guide_legend(ncol = 2), size = guide_legend(ncol = 2)) +
	coord_sf(default_crs = st_crs(4326)) +
	scale_x_continuous(breaks = extended_breaks(10)) +
	theme(legend.position = "bottom", legend.title.position = "top") +
	xlab("") +
	ylab("")

# Plotting
ggarrange(H3N2RePlot, asymNetwork, labels = "AUTO", align = "hv", common.legend = T, legend.grob = get_legend(asymNetwork), legend = "bottom")
ggsave("Figures/PhylodynamicsResults_smaller.svg", width = 9, height = 6) # NOTE: LARGE FILE

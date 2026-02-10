library(FactoMineR)
library(factoextra)
library(cluster)
library(colorspace)
library(lubridate)
library(pbapply)
library(ggrepel)
library(scales)
library(grid)
library(Biostrings)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggnewscale)
library(ggforce)
library(xtable)
library(gtools)
library(svglite)
library(cowplot)
library(phangorn)
library(pheatmap)
library(emmeans)
library(brms)
library(cmdstanr)

colourProvince <- c("Canada" = "#FFE0E2", "BC" = "#FEE1D9", "SK" = "#F1E3F2", "MB" = "#DBE9F5", "ON" = "#D1EFD4", "NB" = "#F8E5C3", "NS" = "#CFEDED")
colourProvinceLine <- c("Canada" = "#D81A21", "British Columbia" = "#C64A1C", "Saskatchewan" = "#92278F", "Manitoba" = "#0369AC",
			"Ontario" = "#2B8737", "New Brunswick" = "#8A600D", "Nova Scotia" = "#367A76")

##### Functions #####
read.aa <- function(file = file.choose(), bin = TRUE, filt = NULL){
	# from https://gist.github.com/shaunpwilkinson/2c4ded3c99a3fe8a08d01c8352bac012
	# Has been modified to create a dataframe if binary not requested.
	# Now also filters the positions for you if requested
  	x <- readLines(file)
  	namelines <- grepl("^>", x)
    	f <- cumsum(namelines)
    	res <- split(x, f)
   	resnames <- sapply(res, function(s) s[1])
      	resnames <- gsub("^>", "", resnames)
        resnames <- gsub("\\|.+", "", resnames) 
        res <- lapply(res, function(s) paste0(s[-1], collapse = ""))
	if(length(filt) > 1){
		res <- lapply(res, function(x) {
			resStrings <- strsplit(x,"") |> unlist()
			resStrings <- resStrings[filt] |> paste0(collapse = "")
			return(resStrings)})
	}
	names(res) <- resnames
	if(bin){
		res <- lapply(res, charToRaw)
		class(res) <- "AAbin"
		return(res)
	}

	res <- lapply(res, function(x) strsplit(x, "") |> unlist()) |> bind_rows() |> t() |> as.data.frame()

	# Finally, this is going to make the parsing easier
	colnames(res) <- paste(gsub(rownames(res)[1],pattern = ".*_", replacement = ""), 1:ncol(res), sep = ":")
	res$Genome <- gsub(":.*|(?<=Reference)_.*", "", rownames(res), perl = T)

	# Now to convert this into a dataframe
	return(res)
}

FluAntigenAnalysis <- function(aaFile, antigenSites, nextStrainData, remove = "none"){

	# The old way made some very key mistakes. This is a fix....
	alignedAA <- read.aa(aaFile, bin = F)
	colnames(alignedAA)[1:(ncol(alignedAA) - 1)] <- 1:(ncol(alignedAA) - 1)
	
	# Removing constant sites so that it is akin to a VCF file
	rownames(alignedAA) <- alignedAA$Genome
	alignedAA <- alignedAA[,-which(colnames(alignedAA) == "Genome")]
	alignedAA <- alignedAA |> t() |> as.data.frame()

	# Switching the / and - to .
	colnames(alignedAA) <- gsub("\\W",".", colnames(alignedAA))

	# Reading in the VCF File
	#vcf <- vcfReading(vcfFile)
	#rownames(vcf) <- vcf$POS
	if(!("none" %in% remove)){
		alignedAA <- alignedAA |> select(-contains(eval(remove)))
	}
	alignedAA <- alignedAA[,sapply(alignedAA, function(x) length(unique(x)) > 1)]

	# Filtering out the positions which are at least 10% gaps/unknown
	noGaps <- apply(alignedAA, MARGIN = 1, function(x){
		      y <- table(x) |> unlist()
		      if(sum(c("X","AB","*") %in% names(y))){
				y <- y/sum(y)
		      		ind <- which(names(y) %in% c("X","AB","*"))
				ifelse(sum(y[ind]) >= 0.1,F,T)
		      }else{
			      return(T)
		      }
	
						  }) 
	alignedAA <- alignedAA[noGaps,]
	
	# If we see that there's any gap, we'll put in the allele for the vaccine
#	noMissing <- pbsapply(alignedAA, cl = 10, function(x){ 
#		      y <- table(x) |> unlist() |> names()
#			ifelse(sum(c("X","AB","*") %in% y), F, T)
#						  })
#	if(any(!noMissing)){
#		colind <- which(!noMissing)
#		if(length(colind) == 1){
#			rowind <- which(alignedAA[,colind] %in% c("X","AB","*"))
#			alignedAA[rowind,colind] <- alignedAA[rowind,1]
#		}else{
#			for(i in colind){
#				rowind <- which(alignedAA[,i] %in% c("X","AB","*"))
#				alignedAA[rowind,i] <- alignedAA[rowind,1]
#			}
#		}
#	
#	}

	# Now we're going to figure out the allele distribution per site
	countPerSite <- alignedAA |> mutate(position = as.numeric(rownames(alignedAA))) |>
		pivot_longer(-position, names_to = "Genome", values_to = "AA") |> 
		left_join(antigenSites, relationship = "many-to-many") |> mutate(Antigen = replace(Antigen, is.na(Antigen), "Other")) |>
		left_join(nextStrainData, relationship = "many-to-many") 

	if("short.clade" %in% colnames(countPerSite)){
		countPerSite <- countPerSite |> unite(CladeDesignation, short.clade, subclade) 
	}else{
		countPerSite <- countPerSite |> rename(CladeDesignation = clade)
	}
	       
	countPerSite <- countPerSite |>  
		mutate(position = as.numeric(position)) |>
	       	group_by(CladeDesignation, Antigen, position) |> count(AA, name = "Genomes") |> 
		pivot_wider(c(Antigen, position, AA), names_from = CladeDesignation, values_from = Genomes, values_fill = 0) |> arrange(Antigen, position, AA) |>
		filter(Antigen != "Other")
	
	# Doing the same thing for the glycosylation sites
	glyco <- nextStrainData |> filter(Genome %in% colnames(alignedAA)) |> mutate(glycosylation = gsub("HA\\d:", "", glycosylation, perl = T))
	
	if("short.clade" %in% colnames(glyco)){
		glyco <- glyco |> unite(CladeDesignation, short.clade, subclade) 
	}else{
		glyco <- glyco |> rename(CladeDesignation = clade)
	}

	glyco <- glyco |> group_by(CladeDesignation) |> count(glycosylation, name = "Genomes")
	
	tmp <- glyco$glycosylation |> sapply(function(x){
					    x <- strsplit(x, ";") |> unlist()
					    site <- gsub(".*:", "", x)
					    names(site) <- gsub(":N.(S|T)", "", x)
					    site <- site[mixedsort(names(site))]
					    return(site)
						  }) |> bind_rows()
	glyco <- glyco |> bind_cols(tmp) |> select(-glycosylation) |> arrange(CladeDesignation, -Genomes)

	return(list("Vars" = alignedAA, "antigenDist" = countPerSite, "GlycoSiteDist" = glyco))
}

dateConverterOntario <- function(x){
	dateSplit <- strsplit(x, " ") |> unlist()
	dateSplit[1] <- months[dateSplit[1]]
	return(paste(dateSplit[c(3,1,2)], sep = "-", collapse = "-"))
}

calculateEllipses <- function(dat){
	# Comes from https://stackoverflow.com/questions/30825387/getting-the-parameters-of-a-data-ellipse-produced-by-the-car-package-in-r
	ell.info <- cov.wt(dat) # Getting the covariance matrix
	eigen.info <- eigen(ell.info$cov) # Need to calculate the eigen vectors
	lengths <- sqrt(eigen.info$values * 2 * qf(0.95,2, nrow(dat) - 1)) # How long are the eigenvectors for a 95%

	if(is.nan(lengths[2])){
		warning("Not enough variablity for an ellipse. Preparing a line instead.\n b will be equal to 1e-6")
		lengths[2] <- 1e-6
	}
	
	# Angle calculation for the major axis
	angle <- atan2(eigen.info$vectors[1,2],eigen.info$vectors[2,2])

	angle = -angle

	df <- data.frame("Centre_x" = ell.info$center[1], "Centre_y" = ell.info$center[2], "a" = lengths[1], "b" = lengths[2], "angle" = angle)
	#df <- data.frame("Centre" = ell.info$center, "xmax" = xmax, "xmin" = xmin, "ymax" = ymax, "ymin" = ymin)
	return(df)
}

pEpitopeComparison <- function(varDataframe, keySitesList, cpus){
	samps <- nrow(varDataframe)

	returnList <- pblapply(1:samps,cl = floor(cpus/3), function(i){
		j = c(i+1):samps

		op <- pboptions(type = "none")
		tmp <- pblapply(j, cl = 3, function(j){
				comps <- lapply(keySitesList, function(x){
				       pEpitope <- sum(varDataframe[i,as.character(x)] != varDataframe[j,as.character(x)]) / length(x)
				       return(data.frame(Start = rownames(varDataframe)[i],End = rownames(varDataframe)[j], Dist = pEpitope))
	 			}) |> bind_rows() |> mutate(AntigenSite = names(keySitesList))
				return(comps)

	 	})
		return(tmp)	
	 })

	return(returnList)
}

plottingLegendHide <- theme(legend.position = "none") # Making it easy to hide the legend when needed

# Setting up some conversions here
months <- c(1:12)
names(months) <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
# Setting up the colours
ontarioColours <- c("#F0454B","#F03093","#F15A22","#39B54A","#3193CC")
darkOntarioColours <- c("#D81A21", "#C64A1C", "#92278F", "#0369AC", "#2B8737")
lightOntarioColours <- c("#FFE0E2", "#FEE1D9", "#F1E3F2", "#DBE9F5", "#D1EFD4")
H1Colours <- c("#DBE9F5","#3193CC","#0369AC", "#FFE0E2","#F0454B", "#D81A21")
colourRamp <- c("#39B54A","#3193CC","#B975B7","#FCAF17","#F0454B","#F03093") 

# This is the shape legend to make things easier
vaccineTypeLegend <- legendGrob(c("Cell", "Egg"), hgap = 0.75, vgap = 0.75, 
                           pch = c(23,22),
			   gp = gpar(col = "black", fill = c("white"),  cex = 2/3), ncol = 1)

vaccineLabelTrans <- list("22.23" = "I", "24" = "II", "25.26" = "III", "26" = "IV")
##### Metadata Info ####
dataSeq <- read.delim("../Nextclade/GisaidNextclade.tab")|> select(seqName,clade, short.clade,subclade,glycosylation) |> rename(Genome = seqName) |> 
	mutate(Genome = gsub("\\|.*","", Genome), Genome = gsub("\\W",".", Genome))

# Want to make sure that the vaccine strains stand out
vaccineIndex <- which(grepl("Egg|Cell", dataSeq$Genome))
vaccineDF <- dataSeq[vaccineIndex,] |>
	separate(Genome, c("TMP", "Vaccine"), remove = F, sep = "\\.", extra = "merge") |> mutate(Vaccine = gsub("_.*", "", Vaccine)) |>
	mutate(short.clade = "Vaccine",clade = vaccineLabelTrans[Vaccine] |> unlist(), subclade = vaccineLabelTrans[Vaccine] |>
	       unlist()) |> select(-Vaccine, -TMP)
dataSeq <- dataSeq[-vaccineIndex,] |> bind_rows(vaccineDF)

# Reading in the location metadata
metadata <- read.delim("../Beast/All_meta_v2.tsv", header = T, col.names = c("Genome", "CollectionDate", "Location")) |> as_tibble() |>
	mutate(Location = replace(Location, Location == "Toronto", "Ontario"), Location = gsub("Nova", "Nova_", Location), Location = gsub("NewB", "New_B", Location))

#### Now to do the H3 using the same method ####
# Let's pull out the glycosylation sites and merge them with the antigenic sites
h3ImportSites <- read.csv("../Data/H3ImportantSitesV3.csv") |> as_tibble() |> rename(position = Position) 
	#mutate(position = as.numeric(gsub("[^0-9]","", position))) |> select(position,Key.Sites,impact.of.change)

detectedGlycoSites <- dataSeq$glycosylation |> strsplit(split = c(":")) |> unlist() |> as.numeric() |> unique()
glycoSites <- data.frame(position = detectedGlycoSites[-1], Antigen = "NGS")
h3ImportSites <- h3ImportSites# |> bind_rows(glycoSites)

h3Results <- FluAntigenAnalysis("../Beast/CanadianSequences/CanadaOnly_Amino.faa",h3ImportSites, dataSeq)

noVars <- h3Results$antigenDist |> count(position, name = "Variants")  |> filter(Variants == 1) |> pull(position) |> unique()
write.csv(h3Results$antigenDist |> filter(!(position %in% noVars)), file = "../Results/H3_SiteMuts.csv", row.names =F, quote = F)
write.csv(h3Results$GlycoSiteDist, file = "../Results/H3_Glyco.csv", quote = F, row.names = F)

# The rest doesn't make the same assumption as above
dataSeq <- read.delim("../Nextclade/GisaidNextclade.tab")|> select(seqName,clade, short.clade,subclade,glycosylation) |> rename(Genome = seqName) |> 
	mutate(Genome = gsub("\\|.*","", Genome), Genome = gsub("\\W",".", Genome))

aaBinaryOnlyAntigen <- read.aa(file = "../Beast/CanadianSequences/CanadaOnly_Amino.faa", bin = T, filt = h3ImportSites$position |> unique())
names(aaBinaryOnlyAntigen) <- gsub("/|-", ".", names(aaBinaryOnlyAntigen))

h3Dist <- dist.ml(aaBinaryOnlyAntigen, model = "FLU") # Calculating the AA distance

fit <- cmdscale(h3Dist, eig = T, k = 4, add = T)
coordH3 <- fit$points |> as_tibble()
coordH3$Genome <- rownames(fit$points)
coordH3 <- coordH3 |> left_join(dataSeq) |> unite(CladeDesignation, short.clade, subclade)
contrib <- round(fit$eig/sum(fit$eig) * 100,2)

# Getting the colours sorted out
cladeCounts <- coordH3 |> filter(!grepl("Egg|Cell", Genome)) |> pull(CladeDesignation) |> sort() |> unique() |> gsub(pattern = "_.*",replacement = "") |> table()
H3Colours <- c(colorRampPalette(colourRamp)(sum(cladeCounts)))
names(H3Colours) <- coordH3 |> filter(!grepl("Egg|Cell", Genome)) |> pull(CladeDesignation) |> sort() |> unique()
H3Colours["2a.3a.1_K"] <- "brown"
names(H3Colours) <- gsub("2a.3a.1_","", names(H3Colours))

H3Legend <- legendGrob(names(H3Colours), hgap = unit(0.5, "lines"), vgap = unit(0.5, "lines"), 
                   pch = 21, gp = gpar(col = "black",fill = H3Colours, pch = 1, cex = 2/3), nrow = 2)

# Getting the vaccines sorted
vaccineIndex <- which(grepl("Egg|Cell", coordH3$Genome))

vaccineDF <- coordH3[vaccineIndex,] |> separate(Genome, c("Genome", "Vaccine"), sep = "\\.", extra = "merge") |> mutate(Vaccine = gsub("_.*", "", Vaccine))
#vaccineColours <- colorRampPalette(c("white", "black"))(vaccineDF$Vaccine |> unique() |> length())
vaccineColours <- c("#ffe119", "#4363d8", "#911eb4", "#46f0f0")
names(vaccineColours) <- vaccineDF$Vaccine |> unique() |> sort()

legendColoursvaccine <- vaccineColours[vaccineDF$Vaccine]
legendShapevaccine <- ifelse(grepl("Egg", coordH3$Genome[vaccineIndex]), 22,23)

# Getting the vaccine colours ready
vaccineLegend <- legendGrob(c("I", "II", "III", "IV"), hgap = 0.75, vgap = 0.75, 
                           pch = 23,
			   gp = gpar(col = "black", fill = vaccineColours,  cex = 2/3), ncol = 2)
## Plotting
savedLegendCustom <- guides(custom = guide_custom(H3Legend,title = "H3", order = 1),
       custom = guide_custom(vaccineLegend, title = "Vaccine Strain", order = 3),# custom = guide_custom(shapeLegend, title = "Type", order = 4), 
       shape = guide_custom(vaccineTypeLegend, title = "Vaccine Type", order = 4),
       colour = guide_none(), fill = guide_none()) 

H3DistAAWAGFilt <- coordH3[-vaccineIndex,] |> ggplot(aes(x = V1, y = V2)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point(data = coordH3[vaccineIndex,], colour = "black", fill = legendColoursvaccine, pch = legendShapevaccine, size = 4, show.legend = F) +
	geom_point(pch = 21, aes(fill = clade), alpha = 0.5, size = 2) + 
	scale_fill_manual(values = H3Colours) +
	theme_classic() +
	xlab(paste0("PCoA 1 (", contrib[1],"%)")) +
	ylab(paste0("PCoA 2 (", contrib[2],"%)")) +
	theme(legend.position = "bottom", legend.title.position = "top") +
	savedLegendCustom

H3DistAAWAGDim34Filt <- coordH3[-vaccineIndex,] |> ggplot(aes(x = V3, y = V4)) +
	geom_hline(yintercept=0, lty = 2, colour = "grey90") + 
	geom_vline(xintercept=0, lty = 2, colour = "grey90") +
	geom_point(data = coordH3[vaccineIndex,], colour = "black", fill = legendColoursvaccine, pch = legendShapevaccine, size = 4, show.legend = F) +
	geom_point(pch = 21, aes(fill = clade)) + 
	scale_fill_manual(values = H3Colours) +
	theme_classic() +
	xlab(paste0("PCoA 3 (", contrib[3],"%)")) +
	ylab(paste0("PCoA 4 (", contrib[4],"%)")) +
	theme(legend.position = "bottom", legend.title.position = "top") 

#segmentsPlotH3 <- ggarrange(H3DistAAWAGFilt + savedLegendCustom, H3DistAAWAGDim34Filt, nrow = 1, align = "hv", labels = "AUTO", common.legend = T, legend = "bottom")

#segmentsPlotH3
#ggsave("../Figures/PCoAAntigensFLU.png", width = 9, height = 6)

# Making a heatmap
legendHeatmap <- coordH3 |> select(Genome, clade)
legendHeatmap <- legendHeatmap[,2] |> as.data.frame()
row.names(legendHeatmap) <- coordH3 |> pull(Genome)
legendHeatmap[vaccineIndex,] <- vaccineDF$Vaccine
coloursList <- list(clade = c(H3Colours, vaccineColours))

#pheatmap(h3Dist, clustering_method = "ward.D2", annotation_row = legendHeatmap, annotation_col = legendHeatmap,
#	 annotation_names_row = F, annotation_names_col = F, annotation_colors = coloursList,
#	 filename = "../Figures/H3AntigenHeatmap.png", width = 6, height = 6, legend_breaks = c(seq(0,0.1,0.02)), legend_labels = c(seq(0,0.1,0.02)))
#dev.off()

h3DistMatrix <- h3Dist |> as.matrix()

cladeList <- split(x = legendHeatmap$clade, f = rownames(legendHeatmap))
#cladesPositionList <- coordH3 |> pull(Genome) |> split(f = legendHeatmap |> as.vector())

idx <- which(lower.tri(h3DistMatrix), arr.ind = T)

AAdistPairwise <- data.frame(Start = colnames(h3DistMatrix)[idx[,1]], End = colnames(h3DistMatrix)[idx[,2]], value = h3DistMatrix[lower.tri(h3DistMatrix)]) |>
	mutate(StartClade = cladeList[Start] |> unlist(), EndClade = cladeList[End] |> unlist())

reformattedPairwise <- AAdistPairwise |> filter(grepl("22.23|24|25.26|26", EndClade)| grepl("22.23|24|25.26|26", StartClade)) |>
	apply(MARGIN = 1, function(x){
	       if(grepl("22.23|24|25.26|26", x[5])){ # Want to put all the vaccines in the same spot
		       newStart <- x[5]
		       newEnd <- x[4]
		       newStartGenome <- x[2]
		       newEndGenome <- x[1]

			x[c(1,2,4,5)] <- c(newStartGenome, newEndGenome, newStart,newEnd)

		       return(x)
			
	       }else{
		       return(x)
	       }
	 }) |> t() |> as.data.frame()
colnames(reformattedPairwise)  <- c("Start", "End", "value", "StartClade", "EndClade")

brmsWAGModel <- reformattedPairwise |> distinct() |> mutate(value = as.numeric(value)) |>
	filter(!grepl("22.23|24|25.26|26", EndClade)) |>
	brm(formula = value+1e-6 ~ 0 + EndClade * StartClade, family = lognormal(), cores = 4, iter = 5e3, backend = "cmdstanr")  # We have a a couple cases of absolute matches. Adding 1e-6 lets us still use the lognormal distribution
estimatedDistances <- brmsWAGModel |> emmeans( ~ StartClade | EndClade, regrid = "response")
estimatedPairs <- estimatedDistances |> pairs()

H3vaccineDistanceWAG <- as.data.frame(estimatedDistances) |>
	mutate(StartClade = vaccineLabelTrans[StartClade] |> unlist(), EndClade = gsub("2a.3a.1_","",EndClade)) |>
	ggplot(aes(y = exp(emmean), ymin = exp(lower.HPD), ymax = exp(upper.HPD), fill = EndClade, x = StartClade, height = 0.5)) +
	geom_col() +
	geom_errorbar(width = 0.5) +
	theme_classic() +
	#ylim(c(0,0.1)) +
	scale_fill_manual(values = coloursList$clade, "Clade") +
	facet_grid(. ~ EndClade) +
	xlab("Vaccine Strain") +
	ylab("Mean FLU Distance") +
	theme(legend.position = "bottom")

H3vaccinePairCompWAG <- as.data.frame(estimatedPairs) |> tibble() |> mutate(contrast = gsub("StartClade","", contrast)) |> #separate(contrast, c("Start", "End"), sep = " - ") |>
	ggplot(aes(x = estimate, xmin = lower.HPD, xmax = upper.HPD, fill = EndClade, y = contrast, height = 0.5)) +
	geom_rect() +
	geom_point() +
	geom_vline(xintercept = 0, lty = 2) +
	theme_classic() +
	scale_fill_manual(values = coloursList$clade, "Clade") +
	facet_grid(. ~ EndClade) +
	ylab("Pairwise Vaccine Comparisons") +
	scale_y_discrete(limits = rev) +
	xlab("Difference in Regresssion Estimates") +
	theme(legend.position = "bottom")

H3vaccinePlotWAG <- ggarrange(H3vaccineDistanceWAG, H3vaccinePairCompWAG, common.legend = T, legend = "bottom", align = "hv", labels = "AUTO", ncol = 1)

### Now to do the provincial tests within only the K Clade
brmsWAGProvincialModel <- reformattedPairwise |> distinct() |> mutate(value = as.numeric(value)) |>
	filter(grepl("K", EndClade)) |> left_join(metadata, by = c("End" = "Genome")) |>
	brm(formula = value+1e-6 ~ 0 + Location * StartClade, family = lognormal(), cores = 4, iter = 5e3, backend = "cmdstanr")  # We have a a couple cases of absolute matches. Adding 1e-6 lets us still use the lognormal distribution

estimatedDistances <- brmsWAGProvincialModel |> emmeans( ~ StartClade | Location, regrid = "response")
estimatedPairs <- estimatedDistances |> pairs()

H3vaccineDistanceWAGProvincial <- as.data.frame(estimatedDistances) |>
	mutate(StartClade = vaccineLabelTrans[StartClade] |> unlist()) |>
	ggplot(aes(y = exp(emmean), ymin = exp(lower.HPD), ymax = exp(upper.HPD), fill = Location, x = StartClade, height = 0.5)) +
	geom_col() +
	geom_errorbar(width = 0.5) +
	theme_classic() +
	#ylim(c(0,0.1)) +
	#scale_fill_manual(values = coloursList$clade, "Clade") +
	facet_grid(. ~ Location) +
	xlab("Vaccine Strain") +
	ylab("Mean FLU Distance") +
	theme(legend.position = "bottom")

#ggsave("../Figures/H3VaccineDistance.png", H3vaccinePlotWAG, width= 9, height = 6)

##### Now to do the same emmeans models using pEpitopes ####
h3AntigenVars <- h3Results$Vars[h3ImportSites$position,] |> t() |> as.data.frame()

antigens <- h3ImportSites$Antigen |> unique()

antigensList <- sapply(antigens, function(x){h3ImportSites |> filter(grepl(x, Antigen)) |> pull(position) })

# Calculating the pEpitope distances
pEpitopeComps <- pEpitopeComparison(h3AntigenVars, antigensList, 16) |> bind_rows() |>
	filter(!is.na(Dist)) |> mutate(StartClade = cladeList[Start] |> unlist(), EndClade = cladeList[End] |> unlist())

reformattedPairwise <- pEpitopeComps |> filter(grepl("22.23|24|25.26|26", EndClade)| grepl("22.23|24|25.26|26", StartClade)) |>
	apply(MARGIN = 1, function(x){
	       if(grepl("22.23|24|25.26|26", x[6])){ # Want to put all the vaccines in the same spot
		       newStart <- x[6]
		       newEnd <- x[5]
		       newStartGenome <- x[2]
		       newEndGenome <- x[1]

			x[c(1,2,5,6)] <- c(newStartGenome, newEndGenome, newStart,newEnd)

		       return(x)
			
	       }else{
		       return(x)
	       }
	 }) |> t() |> as.data.frame()

antigenPairwise <- reformattedPairwise |> filter(grepl("22.23|24|25.26|26", StartClade),!grepl("22.23|24|25.26|26", EndClade)) |>
	mutate(Dist = as.numeric(Dist)) |> filter(AntigenSite != "RBS")

bayesFormula <- bf(Dist + 1e-6 ~ 0 + EndClade * StartClade * AntigenSite, family = Beta())

brmspEpiAntigenModel <- brm(data = antigenPairwise, formula = bayesFormula, cores = 4, iter = 5e3,# Want to use a Beta prior here as pEpitope is limited to [0,1]
	     backend = "cmdstanr", normalize = F) # Trying to speed things up here. Using cmdstanR and removing normalization 

estimatedDistance <- brmspEpiAntigenModel |> emmeans( ~ StartClade | AntigenSite * EndClade, regrid = "response")
brmspEpiAntigenModel |> pairs()

H3vaccineDistanceEpiIndivSites <- as.data.frame(estimatedDistance) |> tibble() |>
	mutate(StartClade = vaccineLabelTrans[StartClade] |> unlist(), EndClade = gsub("2a.3a.1_","",EndClade)) |>
	ggplot(aes(y = response, ymin = lower.HPD, ymax = upper.HPD, fill = EndClade, x = StartClade, height = 0.5)) +
	geom_col() +
	geom_errorbar(width = 0.5) +
	theme_classic() +
	geom_hline(yintercept = 0.2, lty = 2, colour = "black") +
	#ylim(c(0,0.08)) +
	scale_fill_manual(values = coloursList$clade, "Clade") +
	facet_grid(AntigenSite ~ EndClade) +
	xlab("Vaccine Strain") +
	ylab("Mean pEpitope Distance") +
	theme(legend.position = "bottom")

ggsave("../Figures/pEpitopeAllSites.png",H3vaccineDistanceEpiIndivSites, width = 12, height = 9)

# Now to look at the provincial level
bayesFormula <- bf(Dist + 1e-6 ~ 0 + Location * StartClade * AntigenSite, family = Beta())

brmsEpiProvincialModel <- antigenPairwise |> distinct() |> 
	filter(grepl("K", EndClade)) |> left_join(metadata, by = c("End" = "Genome")) |>
	brm(formula = bayesFormula, cores = 4, iter = 5e3, backend = "cmdstanr")  # We have a a couple cases of absolute matches. Adding 1e-6 lets us still use the lognormal distribution

estimatedDistance <- brmsEpiProvincialModel |> emmeans( ~ StartClade | AntigenSite * Location, regrid = "response")
estimatedPairs <- estimatedDistance |> pairs()

H3vaccineDistanceEpiProvincialIndivSites <- as.data.frame(estimatedDistance) |> tibble() |>
	mutate(StartClade = vaccineLabelTrans[StartClade] |> unlist()) |>
	mutate(Location = gsub("_{1,2}"," ", Location, perl = T)) |>
	ggplot(aes(y = response, ymin = lower.HPD, ymax = upper.HPD, fill = Location, x = StartClade, height = 0.5)) +
	geom_col() +
	geom_errorbar(width = 0.5) +
	theme_classic() +
	geom_hline(yintercept = 0.2, lty = 2, colour = "black") +
	#ylim(c(0,0.08)) +
	scale_fill_manual(values = colourProvinceLine, "Province") +
	facet_grid(AntigenSite ~ Location) +
	xlab("Vaccine Strain") +
	ylab("Mean pEpitope Distance") +
	theme(legend.position = "bottom", legend.title.position = "top")

ggsave(filename = "../Figures/KCladePepitopeAll.svg", H3vaccineDistanceEpiProvincialIndivSites, width = 9, height = 6)

# Now to do the same thing, but, we're only selecting the maximum for each group
globalReformattedPairwise <- reformattedPairwise |> group_by(Start, End) |> filter(Dist == max(Dist)) |>
	filter(grepl("22.23|24|25.26|26", StartClade),!grepl("22.23|24|25.26|26", EndClade)) |>
	 mutate(Dist = as.numeric(Dist)) 

bayesFormula <- bf(Dist + 1e-6 ~ 0 + EndClade * StartClade, family = Beta())

brmspEpiGlobalModel <- brm(data = globalReformattedPairwise, formula = bayesFormula, cores = 4, iter = 5e3,# Want to use a Beta prior here as pEpitope is limited to [0,1]
	     backend = "cmdstanr", normalize = F) # Trying to speed things up here. Using cmdstanR and removing normalization 

estimatedDistance <- brmspEpiGlobalModel |> emmeans( ~ StartClade | EndClade, regrid = "response") 
estimatedPairs <- estimatedDistance |> pairs()

H3vaccineDistanceEpi <- as.data.frame(estimatedDistance) |>
	mutate(StartClade = vaccineLabelTrans[StartClade] |> unlist(), EndClade = gsub("2a.3a.1_","",EndClade)) |>
	ggplot(aes(y = response, ymin = lower.HPD, ymax = upper.HPD, fill = EndClade, x = StartClade, height = 0.5)) +
	geom_col() +
	geom_errorbar(width = 0.5) +
	theme_classic() +
	geom_hline(yintercept = 0.2, lty = 2) +
	#ylim(c(0,0.08)) +
	scale_fill_manual(values = coloursList$clade, "Clade") +
	facet_grid(. ~ EndClade) +
	xlab("Vaccine Strain") +
	ylab("Mean pEpitope Distance") +
	theme(legend.position = "bottom")

### Now to do the provincial tests within only the K Clade
brmsEpiGlobalProvincialModel <- globalReformattedPairwise |> distinct() |> mutate(Dist = as.numeric(Dist)) |>
	filter(grepl("K", EndClade)) |> left_join(metadata, by = c("End" = "Genome")) |>
	brm(formula = Dist+1e-6 ~ 0 + Location * StartClade, family = Beta(), cores = 4, iter = 5e3, backend = "cmdstanr")  # We have a a couple cases of absolute matches. Adding 1e-6 lets us still use the lognormal distribution

estimatedDistances <- brmsEpiGlobalProvincialModel |> emmeans( ~ StartClade | Location, regrid = "response")
estimatedPairs <- estimatedDistances |> pairs()

H3vaccineDistanceEpiProvincial <- as.data.frame(estimatedDistances) |>
	mutate(StartClade = vaccineLabelTrans[StartClade] |> unlist()) |>
	ggplot(aes(y = response, ymin = lower.HPD, ymax = upper.HPD, fill = Location, x = StartClade, height = 0.5)) +
	geom_col() +
	geom_errorbar(width = 0.5) +
	theme_classic() +
	#ylim(c(0,0.1)) +
	#scale_fill_manual(values = coloursList$clade, "Clade") +
	facet_grid(. ~ Location) +
	xlab("Vaccine Strain") +
	ylab("Mean pEpitope Distance") +
	theme(legend.position = "bottom")

#H3vaccinePairCompEpi <- as.data.frame(estimatedPairs) |> tibble() |> mutate(contrast = gsub("StartClade","", contrast)) |> #separate(contrast, c("Start", "End"), sep = " - ") |>
#	mutate(EndClade = gsub("2a.3a.1_","",EndClade)) |>
#	ggplot(aes(x = estimate, xmin = lower.HPD, xmax = upper.HPD, fill = EndClade, y = contrast, height = 0.5)) +
#	geom_rect() +
#	geom_point() +
#	geom_vline(xintercept = 0, lty = 2) +
#	theme_classic() +
#	scale_fill_manual(values = coloursList$CladeDesignation, "Clade") +
#	facet_grid(. ~ EndClade) +
#	ylab("Pairwise Vaccine Comparisons") +
#	scale_y_discrete(limits = rev) +
#	xlab("Difference in pEpitope Distance") +
#	theme(legend.position = "bottom")

#H3vaccinePlotEpi <- ggarrange(H3vaccineDistanceEpi, H3vaccinePairCompEpi, common.legend = T, legend = "bottom", align = "hv", labels = "AUTO", ncol = 1)
#
#ggsave("../Figures/H3VaccineDistanceEpi.png", H3vaccinePlotEpi, width= 9, height = 6)
#
#ggarrange(H3vaccineDistanceWAG + theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank()), H3vaccineDistanceEpi, common.legend = T, legend = "bottom", align = "v", labels = "AUTO", ncol = 1)
#ggsave("../Figures/H3VaccineDistances.png", width = 9, height = 6)
#
##### Let's make some figures here ####
ggarrange(H3DistAAWAGFilt, H3vaccineDistanceEpi, nrow = 1, align = "hv", common.legend = T, legend = "bottom", labels = "AUTO")
ggsave("../Figures/H3VaccineDistancesCanadaOnly.svg", width = 9, height = 6)


ggarrange(H3DistAAWAGFilt,
	  ggarrange(H3vaccineDistanceWAG + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
	  H3vaccineDistanceEpi + theme(axis.title.x = element_blank()), nrow = 2, align = "v", legend = "none", common.legend = T, labels = c("B", "C")),
	  common.legend = T, legend = "bottom", align = "v", ncol = 2, labels = c("A", ""), widths = c(2,3))
ggsave("../Figures/H3VaccineDistancesCanadaOnly.svg", width = 12, height = 9)

ggarrange(H3vaccineDistanceWAGProvincial + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
	  H3vaccineDistanceEpiProvincial + theme(axis.title.x = element_blank()), nrow = 2, align = "hv", legend = "none", common.legend = T, labels = c("AUTO"))

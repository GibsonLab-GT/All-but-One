#! /usr/local/pacerepov1/R/3.4.3/bin/Rscript

# PHASE III Pt.3 (add to csv_plots.sh script)

# Overlay LD for each peak

# !!!!!!!!!!!!! HAVE ONLY WRITTEN CODE FOR UP TO 9 PEAKS !!!!!!!!!!!!!!!!!!!!

args = commandArgs(trailingOnly=TRUE)

# create all color palettes to use

rbPalPk1 <- colorRampPalette(c('snow3','blue','dark blue'))
rbPalPk2 <- colorRampPalette(c('snow3','greenyellow','darkgreen'))
rbPalPk3 <- colorRampPalette(c('snow3','khaki1', 'darkgoldenrod3'))
rbPalPk4 <- colorRampPalette(c('snow3','tomato', 'firebrick4'))
rbPalPk5 <- colorRampPalette(c('snow3','orchid', 'purple4'))
rbPalPk6 <- colorRampPalette(c('snow3','tan1', 'chocolate4'))
rbPalPk7 <- colorRampPalette(c('snow3','pink', 'deeppink4'))
rbPalPk8 <- colorRampPalette(c('snow3','orange', 'darkorange4'))
rbPalPk9 <- colorRampPalette(c('snow3','aquamarine', 'aquamarine4'))
rbPalPk10 <- colorRampPalette(c('snow3','lightslateblue', 'slateblue4'))

if (isTRUE(grepl("1", args[2]))) {
	print("1 peak")
	# reads csv files
	library(readr)
	SNP_list <- read.table(args[4])
	univariate <- read_csv(args[7])
	peak1 <- read_csv(args[6])
	conditioned_all <-read_csv(args[5])

	gene <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 20))
	probe <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 1, 7))

	# assign peak colors
	R1 <- rbPalPk1(10)[as.numeric(cut(peak1$'R2', breaks=10))]

	# make univariate plot
	uni_title = paste('LocusZoom', gene, probe, "\n", "Univariate", sep = " ", collapes = NULL)
	uni_output_name = paste(args[3], "_univariate.pdf", sep = "", collapes = NULL)
	pdf(file = uni_output_name)
	plot(univariate$BP/1000000, -log10(univariate$P), xlab='Position (Mb)', ylab='-logPvalue', main=uni_title, pch = 20, col = R1)
	dev.off()

	# make plot for each peak; 1
	peak_SNP <- SNP_list[1, 1]
	pk1_title = paste('LocusZoom', gene, probe, "\n", "Peak 1 SNP", peak_SNP, sep = " ", collapes = NULL)
	pk1_output_name = paste(args[3], "_conditioned_peak1.pdf", sep = "", collapes = NULL)
	pdf(file = pk1_output_name)
	plot(peak1$BP/1000000, -log10(peak1$P), xlab='Position (Mb)', ylab='-logPvalue', main=pk1_title, pch = 20, col = R1)
	dev.off()

	# make plot for conditioned
	cond_title = paste('LocusZoom', gene, probe, "\n", "Conditioned on All", sep = " ", collapes = NULL)
	cond_output_name = paste(args[3], "_conditioned_all.pdf", sep = "", collapes = NULL)
	pdf(file = cond_output_name)
	plot(conditioned_all$BP/1000000, -log10(conditioned_all$P), xlab='Position (Mb)', ylab='-logPvalue', main=cond_title, pch = 20, col = R1)
	dev.off()

} else if (isTRUE(grepl("2", args[2]))) {
	print("2 peaks")
	# reads csv files
	library(readr)
	SNP_list <- read.table(args[4])
	univariate <- read_csv(args[8])
	peak1 <- read_csv(args[6])
	peak2 <- read_csv(args[7])
	conditioned_all <-read_csv(args[5])

	gene <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 20))
	probe <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 1, 7))

	# assign peak colors
	R1 <- rbPalPk1(10)[as.numeric(cut(peak1$'R2', breaks=10))]
	R2 <- rbPalPk2(10)[as.numeric(cut(peak2$'R2', breaks=10))]

	# make univariate plot
	uni_title = paste('LocusZoom', gene, probe, "\n", "Univariate", sep = " ", collapes = NULL)
	uni_output_name = paste(args[3], "_univariate.pdf", sep = "", collapes = NULL)
	pdf(file = uni_output_name)
	plot(univariate$BP/1000000, -log10(univariate$P), xlab='Position (Mb)', ylab='-logPvalue', main=uni_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2'), R1, R2))
	dev.off()

	# make plot for each peak; 

	# Peak 1
	peak1_SNP <- SNP_list[1, 1]
	pk1_title = paste('LocusZoom', gene, probe, "\n", "Peak 1", "SNP", peak1_SNP, sep = " ", collapes = NULL)
	pk1_output_name = paste(args[3], "_conditioned_peak1", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk1_output_name)
	plot((peak1$"BP")/1000000, -log10(peak1$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk1_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2'), R1, R2))
	dev.off()

	# Peak 2
	peak2_SNP <- SNP_list[2, 1]
	pk2_title = paste('LocusZoom', gene, probe, "\n", "Peak 2", "SNP", peak2_SNP, sep = " ", collapes = NULL)
	pk2_output_name = paste(args[3], "_conditioned_peak2", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk2_output_name)
	plot((peak2$"BP")/1000000, -log10(peak2$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk2_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2'), R1, R2))
	dev.off()

	# make plot for conditioned
	cond_title = paste('LocusZoom', gene, probe, "\n", "Conditioned on All", sep = " ", collapes = NULL)
	cond_output_name = paste(args[3], "_conditioned_all.pdf", sep = "", collapes = NULL)
	pdf(file = cond_output_name)
	plot(conditioned_all$BP/1000000, -log10(conditioned_all$P), xlab='Position (Mb)', ylab='-logPvalue', main=cond_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2'), R1, R2))
	dev.off()

} else if (isTRUE(grepl("3", args[2]))) {
	print("3 peaks")
	# reads csv files
	library(readr)
	SNP_list <- read.table(args[4])
	univariate <- read_csv(args[9])
	peak1 <- read_csv(args[6])
	peak2 <- read_csv(args[7])
	peak3 <- read_csv(args[8])
	conditioned_all <-read_csv(args[5])

	gene <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 20))
	probe <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 1, 7))

	# assign peak colors
	R1 <- rbPalPk1(10)[as.numeric(cut(peak1$'R2', breaks=10))]
	R2 <- rbPalPk2(10)[as.numeric(cut(peak2$'R2', breaks=10))]
	R3 <- rbPalPk3(10)[as.numeric(cut(peak3$'R2', breaks=10))]

	# make univariate plot
	uni_title = paste('LocusZoom', gene, probe, "\n", "Univariate", sep = " ", collapes = NULL)
	uni_output_name = paste(args[3], "_univariate.pdf", sep = "", collapes = NULL)
	pdf(file = uni_output_name)
	plot(univariate$BP/1000000, -log10(univariate$P), xlab='Position (Mb)', ylab='-logPvalue', main=uni_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2'), R2, R3)))
	dev.off()

	# make plot for each peak; 

	# Peak 1
	peak1_SNP <- SNP_list[1, 1]
	pk1_title = paste('LocusZoom', gene, probe, "\n", "Peak 1", "SNP", peak1_SNP, sep = " ", collapes = NULL)
	pk1_output_name = paste(args[3], "_conditioned_peak1", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk1_output_name)
	plot((peak1$"BP")/1000000, -log10(peak1$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk1_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2'), R2, R3)))
	dev.off()

	# Peak 2
	peak2_SNP <- SNP_list[2, 1]
	pk2_title = paste('LocusZoom', gene, probe, "\n", "Peak 2", "SNP", peak2_SNP, sep = " ", collapes = NULL)
	pk2_output_name = paste(args[3], "_conditioned_peak2", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk2_output_name)
	plot((peak2$"BP")/1000000, -log10(peak2$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk2_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2'), R2, R3)))
	dev.off()

	# Peak 3
	peak3_SNP <- SNP_list[3, 1]
	pk3_title = paste('LocusZoom', gene, probe, "\n", "Peak 3", "SNP", peak3_SNP, sep = " ", collapes = NULL)
	pk3_output_name = paste(args[3], "_conditioned_peak3", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk3_output_name)
	plot((peak3$"BP")/1000000, -log10(peak3$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk3_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2'), R2, R3)))
	dev.off()

	# make plot for conditioned
	cond_title = paste('LocusZoom', gene, probe, "\n", "Conditioned on All", sep = " ", collapes = NULL)
	cond_output_name = paste(args[3], "_conditioned_all.pdf", sep = "", collapes = NULL)
	pdf(file = cond_output_name)
	plot(conditioned_all$BP/1000000, -log10(conditioned_all$P), xlab='Position (Mb)', ylab='-logPvalue', main=cond_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2'), R2, R3)))
	dev.off()

} else if (isTRUE(grepl("4", args[2]))) {
	print("4 peaks")
	# reads csv files
	library(readr)
	SNP_list <- read.table(args[4])
	univariate <- read_csv(args[10])
	peak1 <- read_csv(args[6])
	peak2 <- read_csv(args[7])
	peak3 <- read_csv(args[8])
	peak4 <- read_csv(args[9])
	conditioned_all <-read_csv(args[5])

	gene <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 20))
	probe <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 1, 7))

	# assign peak colors
	R1 <- rbPalPk1(10)[as.numeric(cut(peak1$'R2', breaks=10))]
	R2 <- rbPalPk2(10)[as.numeric(cut(peak2$'R2', breaks=10))]
	R3 <- rbPalPk3(10)[as.numeric(cut(peak3$'R2', breaks=10))]
	R4 <- rbPalPk4(10)[as.numeric(cut(peak4$'R2', breaks=10))]

	# make univariate plot
	uni_title = paste('LocusZoom', gene, probe, "\n", "Univariate", sep = " ", collapes = NULL)
	uni_output_name = paste(args[3], "_univariate.pdf", sep = "", collapes = NULL)
	pdf(file = uni_output_name)
	plot(univariate$BP/1000000, -log10(univariate$P), xlab='Position (Mb)', ylab='-logPvalue', main=uni_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2'), R3, R4))))
	dev.off()

	# make plot for each peak; 

	# Peak 1
	peak1_SNP <- SNP_list[1, 1]
	pk1_title = paste('LocusZoom', gene, probe, "\n", "Peak 1", "SNP", peak1_SNP, sep = " ", collapes = NULL)
	pk1_output_name = paste(args[3], "_conditioned_peak1", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk1_output_name)
	plot((peak1$"BP")/1000000, -log10(peak1$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk1_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2'), R1, 
		ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2'), R2, 
		ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2'), R3, R4))))
	dev.off()

	# Peak 2
	peak2_SNP <- SNP_list[2, 1]
	pk2_title = paste('LocusZoom', gene, probe, "\n", "Peak 2", "SNP", peak2_SNP, sep = " ", collapes = NULL)
	pk2_output_name = paste(args[3], "_conditioned_peak2", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk2_output_name)
	plot((peak2$"BP")/1000000, -log10(peak2$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk2_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2'), R1, 
		ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2'), R2, 
		ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2'), R3, R4))))
	dev.off()

	# Peak 3
	peak3_SNP <- SNP_list[3, 1]
	pk3_title = paste('LocusZoom', gene, probe, "\n", "Peak 3", "SNP", peak3_SNP, sep = " ", collapes = NULL)
	pk3_output_name = paste(args[3], "_conditioned_peak3", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk3_output_name)
	plot((peak3$"BP")/1000000, -log10(peak3$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk3_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2'), R1, 
		ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2'), R2, 
		ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2'), R3, R4))))
	dev.off()

	# Peak 4
	peak4_SNP <- SNP_list[4, 1]
	pk4_title = paste('LocusZoom', gene, probe, "\n", "Peak 4", "SNP", peak4_SNP, sep = " ", collapes = NULL)
	pk4_output_name = paste(args[3], "_conditioned_peak4", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk4_output_name)
	plot((peak4$"BP")/1000000, -log10(peak4$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk4_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2'), R1, 
		ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2'), R2, 
		ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2'), R3, R4))))
	dev.off()

	# make plot for conditioned
	cond_title = paste('LocusZoom', gene, probe, "\n", "Conditioned on All", sep = " ", collapes = NULL)
	cond_output_name = paste(args[3], "_conditioned_all.pdf", sep = "", collapes = NULL)
	pdf(file = cond_output_name)
	plot(conditioned_all$BP/1000000, -log10(conditioned_all$P), xlab='Position (Mb)', ylab='-logPvalue', main=cond_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2'), R3, R4))))
	dev.off()

} else if (isTRUE(grepl("5", args[2]))) {
	print("5 peaks")
	# reads csv files
	library(readr)
	SNP_list <- read.table(args[4])
	univariate <- read_csv(args[11])
	peak1 <- read_csv(args[6])
	peak2 <- read_csv(args[7])
	peak3 <- read_csv(args[8])
	peak4 <- read_csv(args[9])
	peak5 <- read_csv(args[10])
	conditioned_all <-read_csv(args[5])

	gene <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 20))
	probe <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 1, 7))

	# assign peak colors
	R1 <- rbPalPk1(10)[as.numeric(cut(peak1$'R2', breaks=10))]
	R2 <- rbPalPk2(10)[as.numeric(cut(peak2$'R2', breaks=10))]
	R3 <- rbPalPk3(10)[as.numeric(cut(peak3$'R2', breaks=10))]
	R4 <- rbPalPk4(10)[as.numeric(cut(peak4$'R2', breaks=10))]
	R5 <- rbPalPk5(10)[as.numeric(cut(peak5$'R2', breaks=10))]

	# make univariate plot
	uni_title = paste('LocusZoom', gene, probe, "\n", "Univariate", sep = " ", collapes = NULL)
	uni_output_name = paste(args[3], "_univariate.pdf", sep = "", collapes = NULL)
	pdf(file = uni_output_name)
	plot(univariate$BP/1000000, -log10(univariate$P), xlab='Position (Mb)', ylab='-logPvalue', main=uni_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2'),	R4, R5)))))
	dev.off()

	# make plot for each peak; 

	# Peak 1
	peak1_SNP <- SNP_list[1, 1]
	pk1_title = paste('LocusZoom', gene, probe, "\n", "Peak 1", "SNP", peak1_SNP, sep = " ", collapes = NULL)
	pk1_output_name = paste(args[3], "_conditioned_peak1", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk1_output_name)
	plot((peak1$"BP")/1000000, -log10(peak1$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk1_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2'),	R4, R5)))))
	dev.off()

	# Peak 2
	peak2_SNP <- SNP_list[2, 1]
	pk2_title = paste('LocusZoom', gene, probe, "\n", "Peak 2", "SNP", peak2_SNP, sep = " ", collapes = NULL)
	pk2_output_name = paste(args[3], "_conditioned_peak2", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk2_output_name)
	plot((peak2$"BP")/1000000, -log10(peak2$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk2_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2'),	R4, R5)))))
	dev.off()

	# Peak 3
	peak3_SNP <- SNP_list[3, 1]
	pk3_title = paste('LocusZoom', gene, probe, "\n", "Peak 3", "SNP", peak3_SNP, sep = " ", collapes = NULL)
	pk3_output_name = paste(args[3], "_conditioned_peak3", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk3_output_name)
	plot((peak3$"BP")/1000000, -log10(peak3$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk3_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2'),	R4, R5)))))
	dev.off()

	# Peak 4
	peak4_SNP <- SNP_list[4, 1]
	pk4_title = paste('LocusZoom', gene, probe, "\n", "Peak 4", "SNP", peak4_SNP, sep = " ", collapes = NULL)
	pk4_output_name = paste(args[3], "_conditioned_peak4", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk4_output_name)
	plot((peak4$"BP")/1000000, -log10(peak4$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk4_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2'),	R4, R5)))))
	dev.off()

	# Peak 5
	peak5_SNP <- SNP_list[5, 1]
	pk5_title = paste('LocusZoom', gene, probe, "\n", "Peak 5", "SNP", peak5_SNP, sep = " ", collapes = NULL)
	pk5_output_name = paste(args[3], "_conditioned_peak5", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk5_output_name)
	plot((peak5$"BP")/1000000, -log10(peak5$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk5_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2'),	R4, R5)))))
	dev.off()

	# make plot for conditioned
	cond_title = paste('LocusZoom', gene, probe, "\n", "Conditioned on All", sep = " ", collapes = NULL)
	cond_output_name = paste(args[3], "_conditioned_all.pdf", sep = "", collapes = NULL)
	pdf(file = cond_output_name)
	plot(conditioned_all$BP/1000000, -log10(conditioned_all$P), xlab='Position (Mb)', ylab='-logPvalue', main=cond_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2'),	R4, R5)))))
	dev.off()

} else if (isTRUE(grepl("6", args[2]))) {
	print("6 peaks")
	# reads csv files
	library(readr)
	SNP_list <- read.table(args[4])
	univariate <- read_csv(args[12])
	peak1 <- read_csv(args[6])
	peak2 <- read_csv(args[7])
	peak3 <- read_csv(args[8])
	peak4 <- read_csv(args[9])
	peak5 <- read_csv(args[10])
	peak6 <- read_csv(args[11])
	conditioned_all <-read_csv(args[5])

	gene <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 20))
	probe <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 1, 7))

	# assign peak colors
	R1 <- rbPalPk1(10)[as.numeric(cut(peak1$'R2', breaks=10))]
	R2 <- rbPalPk2(10)[as.numeric(cut(peak2$'R2', breaks=10))]
	R3 <- rbPalPk3(10)[as.numeric(cut(peak3$'R2', breaks=10))]
	R4 <- rbPalPk4(10)[as.numeric(cut(peak4$'R2', breaks=10))]
	R5 <- rbPalPk5(10)[as.numeric(cut(peak5$'R2', breaks=10))]
	R6 <- rbPalPk6(10)[as.numeric(cut(peak6$'R2', breaks=10))]

	# make univariate plot
	uni_title = paste('LocusZoom', gene, probe, "\n", "Univariate", sep = " ", collapes = NULL)
	uni_output_name = paste(args[3], "_univariate.pdf", sep = "", collapes = NULL)
	pdf(file = uni_output_name)
	plot(univariate$BP/1000000, -log10(univariate$P), xlab='Position (Mb)', ylab='-logPvalue', main=uni_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2'), R5, R6))))))
	dev.off()

	# make plot for each peak; 

	# Peak 1
	peak1_SNP <- SNP_list[1, 1]
	pk1_title = paste('LocusZoom', gene, probe, "\n", "Peak 1", "SNP", peak1_SNP, sep = " ", collapes = NULL)
	pk1_output_name = paste(args[3], "_conditioned_peak1", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk1_output_name)
	plot((peak1$"BP")/1000000, -log10(peak1$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk1_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2'), R5, R6))))))
	dev.off()

	# Peak 2
	peak2_SNP <- SNP_list[2, 1]
	pk2_title = paste('LocusZoom', gene, probe, "\n", "Peak 2", "SNP", peak2_SNP, sep = " ", collapes = NULL)
	pk2_output_name = paste(args[3], "_conditioned_peak2", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk2_output_name)
	plot((peak2$"BP")/1000000, -log10(peak2$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk2_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2'), R5, R6))))))
	dev.off()

	# Peak 3
	peak3_SNP <- SNP_list[3, 1]
	pk3_title = paste('LocusZoom', gene, probe, "\n", "Peak 3", "SNP", peak3_SNP, sep = " ", collapes = NULL)
	pk3_output_name = paste(args[3], "_conditioned_peak3", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk3_output_name)
	plot((peak3$"BP")/1000000, -log10(peak3$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk3_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2'), R5, R6))))))
	dev.off()

	# Peak 4
	peak4_SNP <- SNP_list[4, 1]
	pk4_title = paste('LocusZoom', gene, probe, "\n", "Peak 4", "SNP", peak4_SNP, sep = " ", collapes = NULL)
	pk4_output_name = paste(args[3], "_conditioned_peak4", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk4_output_name)
	plot((peak4$"BP")/1000000, -log10(peak4$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk4_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2'), R5, R6))))))
	dev.off()

	# Peak 5
	peak5_SNP <- SNP_list[5, 1]
	pk5_title = paste('LocusZoom', gene, probe, "\n", "Peak 5", "SNP", peak5_SNP, sep = " ", collapes = NULL)
	pk5_output_name = paste(args[3], "_conditioned_peak5", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk5_output_name)
	plot((peak5$"BP")/1000000, -log10(peak5$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk5_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2'), R5, R6))))))
	dev.off()

	# Peak 6
	peak6_SNP <- SNP_list[6, 1]
	pk6_title = paste('LocusZoom', gene, probe, "\n", "Peak 6", "SNP", peak6_SNP, sep = " ", collapes = NULL)
	pk6_output_name = paste(args[3], "_conditioned_peak6", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk6_output_name)
	plot((peak6$"BP")/1000000, -log10(peak6$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk6_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2'& peak5$'R2' > peak6$'R2'), R5, R6))))))
	dev.off()

	# make plot for conditioned
	cond_title = paste('LocusZoom', gene, probe, "\n", "Conditioned on All", sep = " ", collapes = NULL)
	cond_output_name = paste(args[3], "_conditioned_all.pdf", sep = "", collapes = NULL)
	pdf(file = cond_output_name)
	plot(conditioned_all$BP/1000000, -log10(conditioned_all$P), xlab='Position (Mb)', ylab='-logPvalue', main=cond_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2'& peak5$'R2' > peak6$'R2'), R5, R6))))))
	dev.off()

} else if (isTRUE(grepl("7", args[2]))) {
	print("7 peaks")
	# reads csv files
	library(readr)
	SNP_list <- read.table(args[4])
	univariate <- read_csv(args[13])
	peak1 <- read_csv(args[6])
	peak2 <- read_csv(args[7])
	peak3 <- read_csv(args[8])
	peak4 <- read_csv(args[9])
	peak5 <- read_csv(args[10])
	peak6 <- read_csv(args[11])
	peak7 <- read_csv(args[12])
	conditioned_all <-read_csv(args[5])

	gene <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 20))
	probe <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 1, 7))

	# assign peak colors
	R1 <- rbPalPk1(10)[as.numeric(cut(peak1$'R2', breaks=10))]
	R2 <- rbPalPk2(10)[as.numeric(cut(peak2$'R2', breaks=10))]
	R3 <- rbPalPk3(10)[as.numeric(cut(peak3$'R2', breaks=10))]
	R4 <- rbPalPk4(10)[as.numeric(cut(peak4$'R2', breaks=10))]
	R5 <- rbPalPk5(10)[as.numeric(cut(peak5$'R2', breaks=10))]
	R6 <- rbPalPk6(10)[as.numeric(cut(peak6$'R2', breaks=10))]
	R7 <- rbPalPk7(10)[as.numeric(cut(peak6$'R2', breaks=10))]

	# make univariate plot
	uni_title = paste('LocusZoom', gene, probe, "\n", "Univariate", sep = " ", collapes = NULL)
	uni_output_name = paste(args[3], "_univariate.pdf", sep = "", collapes = NULL)
	pdf(file = uni_output_name)
	plot(univariate$BP/1000000, -log10(univariate$P), xlab='Position (Mb)', ylab='-logPvalue', main=uni_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2'), R6, R7)))))))
	dev.off()

	# make plot for each peak; 

	# Peak 1
	peak1_SNP <- SNP_list[1, 1]
	pk1_title = paste('LocusZoom', gene, probe, "\n", "Peak 1", "SNP", peak1_SNP, sep = " ", collapes = NULL)
	pk1_output_name = paste(args[3], "_conditioned_peak1", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk1_output_name)
	plot((peak1$"BP")/1000000, -log10(peak1$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk1_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2'), R6, R7)))))))
	dev.off()

	# Peak 2
	peak2_SNP <- SNP_list[2, 1]
	pk2_title = paste('LocusZoom', gene, probe, "\n", "Peak 2", "SNP", peak2_SNP, sep = " ", collapes = NULL)
	pk2_output_name = paste(args[3], "_conditioned_peak2", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk2_output_name)
	plot((peak2$"BP")/1000000, -log10(peak2$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk2_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2'), R6, R7)))))))
	dev.off()

	# Peak 3
	peak3_SNP <- SNP_list[3, 1]
	pk3_title = paste('LocusZoom', gene, probe, "\n", "Peak 3", "SNP", peak3_SNP, sep = " ", collapes = NULL)
	pk3_output_name = paste(args[3], "_conditioned_peak3", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk3_output_name)
	plot((peak3$"BP")/1000000, -log10(peak3$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk3_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2'), R6, R7)))))))
	dev.off()

	# Peak 4
	peak4_SNP <- SNP_list[4, 1]
	pk4_title = paste('LocusZoom', gene, probe, "\n", "Peak 4", "SNP", peak4_SNP, sep = " ", collapes = NULL)
	pk4_output_name = paste(args[3], "_conditioned_peak4", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk4_output_name)
	plot((peak4$"BP")/1000000, -log10(peak4$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk4_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2'), R6, R7)))))))
	dev.off()

	# Peak 5
	peak5_SNP <- SNP_list[5, 1]
	pk5_title = paste('LocusZoom', gene, probe, "\n", "Peak 5", "SNP", peak5_SNP, sep = " ", collapes = NULL)
	pk5_output_name = paste(args[3], "_conditioned_peak5", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk5_output_name)
	plot((peak5$"BP")/1000000, -log10(peak5$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk5_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2'), R6, R7)))))))
	dev.off()

	# Peak 6
	peak6_SNP <- SNP_list[6, 1]
	pk6_title = paste('LocusZoom', gene, probe, "\n", "Peak 6", "SNP", peak6_SNP, sep = " ", collapes = NULL)
	pk6_output_name = paste(args[3], "_conditioned_peak6", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk6_output_name)
	plot((peak6$"BP")/1000000, -log10(peak6$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk6_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2'), R6, R7)))))))
	dev.off()

	# Peak 7
	peak7_SNP <- SNP_list[7, 1]
	pk7_title = paste('LocusZoom', gene, probe, "\n", "Peak 7", "SNP", peak7_SNP, sep = " ", collapes = NULL)
	pk7_output_name = paste(args[3], "_conditioned_peak7", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk7_output_name)
	plot((peak7$"BP")/1000000, -log10(peak7$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk7_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2'), R6, R7)))))))
	dev.off()

	# make plot for conditioned
	cond_title = paste('LocusZoom', gene, probe, "\n", "Conditioned on All", sep = " ", collapes = NULL)
	cond_output_name = paste(args[3], "_conditioned_all.pdf", sep = "", collapes = NULL)
	pdf(file = cond_output_name)
	plot(conditioned_all$BP/1000000, -log10(conditioned_all$P), xlab='Position (Mb)', ylab='-logPvalue', main=cond_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2'), R6, R7)))))))
	dev.off()

} else if (isTRUE(grepl("8", args[2]))) {
	print("8 peaks")
	# reads csv files
	library(readr)
	SNP_list <- read.table(args[4])
	univariate <- read_csv(args[14])
	peak1 <- read_csv(args[6])
	peak2 <- read_csv(args[7])
	peak3 <- read_csv(args[8])
	peak4 <- read_csv(args[9])
	peak5 <- read_csv(args[10])
	peak6 <- read_csv(args[11])
	peak7 <- read_csv(args[12])
	peak8 <- read_csv(args[13])
	conditioned_all <-read_csv(args[5])

	gene <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 20))
	probe <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 1, 7))

	# assign peak colors
	R1 <- rbPalPk1(10)[as.numeric(cut(peak1$'R2', breaks=10))]
	R2 <- rbPalPk2(10)[as.numeric(cut(peak2$'R2', breaks=10))]
	R3 <- rbPalPk3(10)[as.numeric(cut(peak3$'R2', breaks=10))]
	R4 <- rbPalPk4(10)[as.numeric(cut(peak4$'R2', breaks=10))]
	R5 <- rbPalPk5(10)[as.numeric(cut(peak5$'R2', breaks=10))]
	R6 <- rbPalPk6(10)[as.numeric(cut(peak6$'R2', breaks=10))]
	R7 <- rbPalPk7(10)[as.numeric(cut(peak7$'R2', breaks=10))]
	R8 <- rbPalPk8(10)[as.numeric(cut(peak8$'R2', breaks=10))]

	# make univariate plot
	uni_title = paste('LocusZoom', gene, probe, "\n", "Univariate", sep = " ", collapes = NULL)
	uni_output_name = paste(args[3], "_univariate.pdf", sep = "", collapes = NULL)
	pdf(file = uni_output_name)
	plot(univariate$BP/1000000, -log10(univariate$P), xlab='Position (Mb)', ylab='-logPvalue', main=uni_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2'), R7, R8))))))))
	dev.off()

	# make plot for each peak; 

	# Peak 1
	peak1_SNP <- SNP_list[1, 1]
	pk1_title = paste('LocusZoom', gene, probe, "\n", "Peak 1", "SNP", peak1_SNP, sep = " ", collapes = NULL)
	pk1_output_name = paste(args[3], "_conditioned_peak1", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk1_output_name)
	plot((peak1$"BP")/1000000, -log10(peak1$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk1_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2'), R7, R8))))))))
	dev.off()

	# Peak 2
	peak2_SNP <- SNP_list[2, 1]
	pk2_title = paste('LocusZoom', gene, probe, "\n", "Peak 2", "SNP", peak2_SNP, sep = " ", collapes = NULL)
	pk2_output_name = paste(args[3], "_conditioned_peak2", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk2_output_name)
	plot((peak2$"BP")/1000000, -log10(peak2$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk2_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2'), R7, R8))))))))
	dev.off()

	# Peak 3
	peak3_SNP <- SNP_list[3, 1]
	pk3_title = paste('LocusZoom', gene, probe, "\n", "Peak 3", "SNP", peak3_SNP, sep = " ", collapes = NULL)
	pk3_output_name = paste(args[3], "_conditioned_peak3", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk3_output_name)
	plot((peak3$"BP")/1000000, -log10(peak3$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk3_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2'), R7, R8))))))))
	dev.off()

	# Peak 4
	peak4_SNP <- SNP_list[4, 1]
	pk4_title = paste('LocusZoom', gene, probe, "\n", "Peak 4", "SNP", peak4_SNP, sep = " ", collapes = NULL)
	pk4_output_name = paste(args[3], "_conditioned_peak4", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk4_output_name)
	plot((peak4$"BP")/1000000, -log10(peak4$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk4_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2'), R7, R8))))))))
	dev.off()

	# Peak 5
	peak5_SNP <- SNP_list[5, 1]
	pk5_title = paste('LocusZoom', gene, probe, "\n", "Peak 5", "SNP", peak5_SNP, sep = " ", collapes = NULL)
	pk5_output_name = paste(args[3], "_conditioned_peak5", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk5_output_name)
	plot((peak5$"BP")/1000000, -log10(peak5$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk5_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2'), R7, R8))))))))
	dev.off()

	# Peak 6
	peak6_SNP <- SNP_list[6, 1]
	pk6_title = paste('LocusZoom', gene, probe, "\n", "Peak 6", "SNP", peak6_SNP, sep = " ", collapes = NULL)
	pk6_output_name = paste(args[3], "_conditioned_peak6", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk6_output_name)
	plot((peak6$"BP")/1000000, -log10(peak6$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk6_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2'), R7, R8))))))))
	dev.off()

	# Peak 7
	peak7_SNP <- SNP_list[7, 1]
	pk7_title = paste('LocusZoom', gene, probe, "\n", "Peak 7", "SNP", peak7_SNP, sep = " ", collapes = NULL)
	pk7_output_name = paste(args[3], "_conditioned_peak7", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk7_output_name)
	plot((peak7$"BP")/1000000, -log10(peak7$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk7_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2'), R7, R8))))))))
	dev.off()

	# Peak 8
	peak8_SNP <- SNP_list[8, 1]
	pk8_title = paste('LocusZoom', gene, probe, "\n", "Peak 8", "SNP", peak8_SNP, sep = " ", collapes = NULL)
	pk8_output_name = paste(args[3], "_conditioned_peak8", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk8_output_name)
	plot((peak8$"BP")/1000000, -log10(peak8$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk8_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2'), R7, R8))))))))
	dev.off()

	# make plot for conditioned
	cond_title = paste('LocusZoom', gene, probe, "\n", "Conditioned on All", sep = " ", collapes = NULL)
	cond_output_name = paste(args[3], "_conditioned_all.pdf", sep = "", collapes = NULL)
	pdf(file = cond_output_name)
	plot(conditioned_all$BP/1000000, -log10(conditioned_all$P), xlab='Position (Mb)', ylab='-logPvalue', main=cond_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2'), R7, R8))))))))
	dev.off()

} else if (isTRUE(grepl("9", args[2]))) {
	print("9 peaks")
	# reads csv files
	library(readr)
	SNP_list <- read.table(args[4])
	univariate <- read_csv(args[15])
	peak1 <- read_csv(args[6])
	peak2 <- read_csv(args[7])
	peak3 <- read_csv(args[8])
	peak4 <- read_csv(args[9])
	peak5 <- read_csv(args[10])
	peak6 <- read_csv(args[11])
	peak7 <- read_csv(args[12])
	peak8 <- read_csv(args[13])
	peak9 <- read_csv(args[14])
	conditioned_all <-read_csv(args[5])

	gene <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 20))
	probe <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 1, 7))

	# assign peak colors
	R1 <- rbPalPk1(10)[as.numeric(cut(peak1$'R2', breaks=10))]
	R2 <- rbPalPk2(10)[as.numeric(cut(peak2$'R2', breaks=10))]
	R3 <- rbPalPk3(10)[as.numeric(cut(peak3$'R2', breaks=10))]
	R4 <- rbPalPk4(10)[as.numeric(cut(peak4$'R2', breaks=10))]
	R5 <- rbPalPk5(10)[as.numeric(cut(peak5$'R2', breaks=10))]
	R6 <- rbPalPk6(10)[as.numeric(cut(peak6$'R2', breaks=10))]
	R7 <- rbPalPk7(10)[as.numeric(cut(peak7$'R2', breaks=10))]
	R8 <- rbPalPk8(10)[as.numeric(cut(peak8$'R2', breaks=10))]
	R9 <- rbPalPk9(10)[as.numeric(cut(peak9$'R2', breaks=10))]

	# make univariate plot
	uni_title = paste('LocusZoom', gene, probe, "\n", "Univariate", sep = " ", collapes = NULL)
	uni_output_name = paste(args[3], "_univariate.pdf", sep = "", collapes = NULL)
	pdf(file = uni_output_name)
	plot(univariate$BP/1000000, -log10(univariate$P), xlab='Position (Mb)', ylab='-logPvalue', main=uni_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

	# make plot for each peak; 

	# Peak 1
	peak1_SNP <- SNP_list[1, 1]
	pk1_title = paste('LocusZoom', gene, probe, "\n", "Peak 1", "SNP", peak1_SNP, sep = " ", collapes = NULL)
	pk1_output_name = paste(args[3], "_conditioned_peak1", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk1_output_name)
	plot((peak1$"BP")/1000000, -log10(peak1$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk1_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

	# Peak 2
	peak2_SNP <- SNP_list[2, 1]
	pk2_title = paste('LocusZoom', gene, probe, "\n", "Peak 2", "SNP", peak2_SNP, sep = " ", collapes = NULL)
	pk2_output_name = paste(args[3], "_conditioned_peak2", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk2_output_name)
	plot((peak2$"BP")/1000000, -log10(peak2$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk2_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

	# Peak 3
	peak3_SNP <- SNP_list[3, 1]
	pk3_title = paste('LocusZoom', gene, probe, "\n", "Peak 3", "SNP", peak3_SNP, sep = " ", collapes = NULL)
	pk3_output_name = paste(args[3], "_conditioned_peak3", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk3_output_name)
	plot((peak3$"BP")/1000000, -log10(peak3$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk3_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

	# Peak 4
	peak4_SNP <- SNP_list[4, 1]
	pk4_title = paste('LocusZoom', gene, probe, "\n", "Peak 4", "SNP", peak4_SNP, sep = " ", collapes = NULL)
	pk4_output_name = paste(args[3], "_conditioned_peak4", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk4_output_name)
	plot((peak4$"BP")/1000000, -log10(peak4$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk4_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

	# Peak 5
	peak5_SNP <- SNP_list[5, 1]
	pk5_title = paste('LocusZoom', gene, probe, "\n", "Peak 5", "SNP", peak5_SNP, sep = " ", collapes = NULL)
	pk5_output_name = paste(args[3], "_conditioned_peak5", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk5_output_name)
	plot((peak5$"BP")/1000000, -log10(peak5$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk5_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

	# Peak 6
	peak6_SNP <- SNP_list[6, 1]
	pk6_title = paste('LocusZoom', gene, probe, "\n", "Peak 6", "SNP", peak6_SNP, sep = " ", collapes = NULL)
	pk6_output_name = paste(args[3], "_conditioned_peak6", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk6_output_name)
	plot((peak6$"BP")/1000000, -log10(peak6$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk6_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

	# Peak 7
	peak7_SNP <- SNP_list[7, 1]
	pk7_title = paste('LocusZoom', gene, probe, "\n", "Peak 7", "SNP", peak7_SNP, sep = " ", collapes = NULL)
	pk7_output_name = paste(args[3], "_conditioned_peak7", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk7_output_name)
	plot((peak7$"BP")/1000000, -log10(peak7$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk7_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

	# Peak 8
	peak8_SNP <- SNP_list[8, 1]
	pk8_title = paste('LocusZoom', gene, probe, "\n", "Peak 8", "SNP", peak8_SNP, sep = " ", collapes = NULL)
	pk8_output_name = paste(args[3], "_conditioned_peak8", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk8_output_name)
	plot((peak8$"BP")/1000000, -log10(peak8$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk8_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

	# Peak 9
	peak9_SNP <- SNP_list[9, 1]
	pk9_title = paste('LocusZoom', gene, probe, "\n", "Peak 9", "SNP", peak9_SNP, sep = " ", collapes = NULL)
	pk9_output_name = paste(args[3], "_conditioned_peak9", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk9_output_name)
	plot((peak9$"BP")/1000000, -log10(peak9$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk9_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

	# make plot for conditioned
	cond_title = paste('LocusZoom', gene, probe, "\n", "Conditioned on All", sep = " ", collapes = NULL)
	cond_output_name = paste(args[3], "_conditioned_all.pdf", sep = "", collapes = NULL)
	pdf(file = cond_output_name)
	plot(conditioned_all$BP/1000000, -log10(conditioned_all$P), xlab='Position (Mb)', ylab='-logPvalue', main=cond_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2'), R8, R9)))))))))
	dev.off()

} else if (isTRUE(grepl("10", args[2]))) {
	print("10 peaks")
	# reads csv files
	library(readr)
	SNP_list <- read.table(args[4])
	univariate <- read_csv(args[16])
	peak1 <- read_csv(args[6])
	peak2 <- read_csv(args[7])
	peak3 <- read_csv(args[8])
	peak4 <- read_csv(args[9])
	peak5 <- read_csv(args[10])
	peak6 <- read_csv(args[11])
	peak7 <- read_csv(args[12])
	peak8 <- read_csv(args[13])
	peak9 <- read_csv(args[14])
	peak10 <- read_csv(args[15])
	conditioned_all <-read_csv(args[5])

	gene <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 20))
	probe <- stringi::stri_reverse(substring(stringi::stri_reverse(args[1]), 1, 7))

	# assign peak colors
	R1 <- rbPalPk1(10)[as.numeric(cut(peak1$'R2', breaks=10))]
	R2 <- rbPalPk2(10)[as.numeric(cut(peak2$'R2', breaks=10))]
	R3 <- rbPalPk3(10)[as.numeric(cut(peak3$'R2', breaks=10))]
	R4 <- rbPalPk4(10)[as.numeric(cut(peak4$'R2', breaks=10))]
	R5 <- rbPalPk5(10)[as.numeric(cut(peak5$'R2', breaks=10))]
	R6 <- rbPalPk6(10)[as.numeric(cut(peak6$'R2', breaks=10))]
	R7 <- rbPalPk7(10)[as.numeric(cut(peak7$'R2', breaks=10))]
	R8 <- rbPalPk8(10)[as.numeric(cut(peak8$'R2', breaks=10))]
	R9 <- rbPalPk9(10)[as.numeric(cut(peak9$'R2', breaks=10))]
	R10 <- rbPalPk10(10)[as.numeric(cut(peak10$'R2', breaks=10))]

	# make univariate plot
	uni_title = paste('LocusZoom', gene, probe, "\n", "Univariate", sep = " ", collapes = NULL)
	uni_output_name = paste(args[3], "_univariate.pdf", sep = "", collapes = NULL)
	pdf(file = uni_output_name)
	plot(univariate$BP/1000000, -log10(univariate$P), xlab='Position (Mb)', ylab='-logPvalue', main=uni_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()

	# make plot for each peak; 

	# Peak 1
	peak1_SNP <- SNP_list[1, 1]
	pk1_title = paste('LocusZoom', gene, probe, "\n", "Peak 1", "SNP", peak1_SNP, sep = " ", collapes = NULL)
	pk1_output_name = paste(args[3], "_conditioned_peak1", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk1_output_name)
	plot((peak1$"BP")/1000000, -log10(peak1$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk1_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()

	# Peak 2
	peak2_SNP <- SNP_list[2, 1]
	pk2_title = paste('LocusZoom', gene, probe, "\n", "Peak 2", "SNP", peak2_SNP, sep = " ", collapes = NULL)
	pk2_output_name = paste(args[3], "_conditioned_peak2", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk2_output_name)
	plot((peak2$"BP")/1000000, -log10(peak2$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk2_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()

	# Peak 3
	peak3_SNP <- SNP_list[3, 1]
	pk3_title = paste('LocusZoom', gene, probe, "\n", "Peak 3", "SNP", peak3_SNP, sep = " ", collapes = NULL)
	pk3_output_name = paste(args[3], "_conditioned_peak3", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk3_output_name)
	plot((peak3$"BP")/1000000, -log10(peak3$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk3_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()

	# Peak 4
	peak4_SNP <- SNP_list[4, 1]
	pk4_title = paste('LocusZoom', gene, probe, "\n", "Peak 4", "SNP", peak4_SNP, sep = " ", collapes = NULL)
	pk4_output_name = paste(args[3], "_conditioned_peak4", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk4_output_name)
	plot((peak4$"BP")/1000000, -log10(peak4$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk4_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()

	# Peak 5
	peak5_SNP <- SNP_list[5, 1]
	pk5_title = paste('LocusZoom', gene, probe, "\n", "Peak 5", "SNP", peak5_SNP, sep = " ", collapes = NULL)
	pk5_output_name = paste(args[3], "_conditioned_peak5", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk5_output_name)
	plot((peak5$"BP")/1000000, -log10(peak5$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk5_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()

	# Peak 6
	peak6_SNP <- SNP_list[6, 1]
	pk6_title = paste('LocusZoom', gene, probe, "\n", "Peak 6", "SNP", peak6_SNP, sep = " ", collapes = NULL)
	pk6_output_name = paste(args[3], "_conditioned_peak6", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk6_output_name)
	plot((peak6$"BP")/1000000, -log10(peak6$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk6_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()

	# Peak 7
	peak7_SNP <- SNP_list[7, 1]
	pk7_title = paste('LocusZoom', gene, probe, "\n", "Peak 7", "SNP", peak7_SNP, sep = " ", collapes = NULL)
	pk7_output_name = paste(args[3], "_conditioned_peak7", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk7_output_name)
	plot((peak7$"BP")/1000000, -log10(peak7$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk7_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()

	# Peak 8
	peak8_SNP <- SNP_list[8, 1]
	pk8_title = paste('LocusZoom', gene, probe, "\n", "Peak 8", "SNP", peak8_SNP, sep = " ", collapes = NULL)
	pk8_output_name = paste(args[3], "_conditioned_peak8", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk8_output_name)
	plot((peak8$"BP")/1000000, -log10(peak8$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk8_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()

	# Peak 9
	peak9_SNP <- SNP_list[9, 1]
	pk9_title = paste('LocusZoom', gene, probe, "\n", "Peak 9", "SNP", peak9_SNP, sep = " ", collapes = NULL)
	pk9_output_name = paste(args[3], "_conditioned_peak9", ".pdf", sep = "", collapes = NULL)
	pdf(file = pk9_output_name)
	plot((peak9$"BP")/1000000, -log10(peak9$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk9_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()

	# Peak 10
	peak10_SNP <- SNP_list[10, 1]
	pk10_title = paste('LocusZoom', gene, probe, "\n", "Peak 10", "SNP", peak10_SNP, sep = " ", collapes = NULL)
	pk10_output_name = paste(args[3], "_conditioned_peak10", ".pdf", sep = "", collapes = NULL)
	pdf(pk10_output_name)
	plot((peak10$"BP")/1000000, -log10(peak10$"P"), xlab='Position (Mb)', ylab='-logPvalue', main=pk10_title, pch = 20,  
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()
	
	# make plot for conditioned
	cond_title = paste('LocusZoom', gene, probe, "\n", "Conditioned on All", sep = " ", collapes = NULL)
	cond_output_name = paste(args[3], "_conditioned_all.pdf", sep = "", collapes = NULL)
	pdf(file = cond_output_name)
	plot(conditioned_all$BP/1000000, -log10(conditioned_all$P), xlab='Position (Mb)', ylab='-logPvalue', main=cond_title, pch = 20, 
		col = ifelse((peak1$'R2' > peak2$'R2' & peak1$'R2' > peak3$'R2' & peak1$'R2' > peak4$'R2' & peak1$'R2' > peak5$'R2' & peak1$'R2' > peak6$'R2' & peak1$'R2' > peak7$'R2' & peak1$'R2' > peak8$'R2' & peak1$'R2' > peak9$'R2' & peak1$'R2' > peak10$'R2'), R1, 
			ifelse((peak2$'R2' > peak1$'R2' & peak2$'R2' > peak3$'R2' & peak2$'R2' > peak4$'R2' & peak2$'R2' > peak5$'R2' & peak2$'R2' > peak6$'R2' & peak2$'R2' > peak7$'R2' & peak2$'R2' > peak8$'R2' & peak2$'R2' > peak9$'R2' & peak2$'R2' > peak10$'R2'), R2, 
			ifelse((peak3$'R2' > peak1$'R2' & peak3$'R2' > peak2$'R2' & peak3$'R2' > peak4$'R2' & peak3$'R2' > peak5$'R2' & peak3$'R2' > peak6$'R2' & peak3$'R2' > peak7$'R2' & peak3$'R2' > peak8$'R2' & peak3$'R2' > peak9$'R2' & peak3$'R2' > peak10$'R2'), R3, 
			ifelse((peak4$'R2' > peak1$'R2' & peak4$'R2' > peak2$'R2' & peak4$'R2' > peak3$'R2' & peak4$'R2' > peak5$'R2' & peak4$'R2' > peak6$'R2' & peak4$'R2' > peak7$'R2' & peak4$'R2' > peak8$'R2' & peak4$'R2' > peak9$'R2' & peak4$'R2' > peak10$'R2'), R4, 
			ifelse((peak5$'R2' > peak1$'R2' & peak5$'R2' > peak2$'R2' & peak5$'R2' > peak3$'R2' & peak5$'R2' > peak4$'R2' & peak5$'R2' > peak6$'R2' & peak5$'R2' > peak7$'R2' & peak5$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak5$'R2' > peak10$'R2'), R5, 
			ifelse((peak6$'R2' > peak1$'R2' & peak6$'R2' > peak2$'R2' & peak6$'R2' > peak3$'R2' & peak6$'R2' > peak4$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak7$'R2' & peak6$'R2' > peak8$'R2' & peak5$'R2' > peak9$'R2' & peak6$'R2' > peak10$'R2'), R6, 
			ifelse((peak7$'R2' > peak1$'R2' & peak7$'R2' > peak2$'R2' & peak7$'R2' > peak3$'R2' & peak7$'R2' > peak4$'R2' & peak7$'R2' > peak5$'R2' & peak7$'R2' > peak6$'R2' & peak7$'R2' > peak8$'R2' & peak7$'R2' > peak9$'R2' & peak7$'R2' > peak10$'R2'), R7, 
			ifelse((peak8$'R2' > peak1$'R2' & peak8$'R2' > peak2$'R2' & peak8$'R2' > peak3$'R2' & peak8$'R2' > peak4$'R2' & peak8$'R2' > peak5$'R2' & peak8$'R2' > peak6$'R2' & peak8$'R2' > peak7$'R2' & peak8$'R2' > peak9$'R2' & peak8$'R2' > peak10$'R2'), R8,
			ifelse((peak9$'R2' > peak1$'R2' & peak9$'R2' > peak2$'R2' & peak9$'R2' > peak3$'R2' & peak9$'R2' > peak4$'R2' & peak9$'R2' > peak5$'R2' & peak9$'R2' > peak6$'R2' & peak9$'R2' > peak7$'R2' & peak9$'R2' > peak8$'R2' & peak9$'R2' > peak10$'R2'), R9, R10))))))))))
	dev.off()
}



 




























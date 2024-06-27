library (ggplot2)
library(dplyr)
library(dada2)
library(phyloseq)
library(vegan)
library(DECIPHER)
library(Biostrings)


main_dir <- dirname(rstudioapi::getSourceEditorContext()$path) 
setwd(main_dir)

path <- './raw_data/Central'
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,140),
                     maxN=0, maxEE=c(2,4), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE)
out

path <- './raw_data/Central/filtered'
list.files(path)
fnFs <- sort(list.files(path, pattern="_F_filt.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R_filt.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
seqtab <- makeSequenceTable(mergers)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
getN <- function(x) sum(getUniques(x))

out <- out[apply(out, 1,function(x) all(x[2]>0)),]
out

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

#taxa <- assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fa", multithread=TRUE)
#taxa <- addSpecies(taxa, "silva_species_assignment_v132.fa")

asv_tab <- t(seqtab.nochim)
write.csv(asv_tab, 'asv_tab_Central.csv', quote=FALSE)

taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138_train_set.fa.gz", multithread=TRUE)

bad_taxa = as.character(rownames(taxa))[(taxa[, 1] == 'Eukaryota' & 
                                           is.na(taxa[, 3])) | is.na(taxa[, 2])]
print(length(bad_taxa))
print(ncol(seqtab.nochim))
seqtab.nochim = seqtab.nochim[, sapply(as.character(colnames(seqtab.nochim)), 
                                       function(x) !(x %in% bad_taxa))]
print(nrow(seqtab.nochim))

#taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138_wSpecies_train_set.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "silva_species_assignment_v138.fa.gz", verbose=TRUE, allowMultiple=T)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
taxa_print_df <- data.frame(taxa.print)

head(taxa_print_df)

write.csv(taxa, "taxa_df_Central.csv", quote=FALSE)
write.csv(taxa_print_df, "taxa_print_Central.csv", row.names=FALSE, quote=FALSE)

colnames(seqtab.nochim) <-  as.character(sapply(colnames(seqtab.nochim), function(x) gsub('NNNNNNNNNN', '', x)))

samples_data <- data.frame(SampleID=sample.names)
rownames(samples_data) <- sample.names

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samples_data), 
               tax_table(taxa))


rich = estimate_richness(ps)
colnames(rich)[0] <- 'Sample'
rich <- rich %>% mutate_if(is.numeric, round, digits = 3)
rich

plot_richness(ps, measures=c("Chao1", "Shannon", "Simpson"))



#ord.pcoa.bray <- ordinate(ps.prop, method="PCoA", distance="bray", na.rm = TRUE)

ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, title="Bray NMDS")

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill="Family")



#abundances_table <- function(ps, seqtab, rank){
#rank_ps <- tax_glom(ps, rank)
#otu_df <- as.data.frame(t(otu_table(rank_ps)))
# taxa_table <- tax_table(rank_ps)
#  taxa_table <- taxa_table[,colSums(is.na(taxa_table))<nrow(taxa_table)] # all rows are saved here
# rownames(otu_df) <- apply(taxa_table, 1, paste, collapse=";")
#if (rank=='Species'){
#   otu_df["taxa"] <- apply(taxa_table[, c('Genus', 'Species')], 1, paste, collapse="_")
#  }
#  else {
#    otu_df["taxa"] <- taxa_table[, rank]
#  }
#  otu_df <- otu_df %>% group_by(taxa) %>% summarise_all(funs(sum))
#  otu_df <- as.data.frame(otu_df)
#  rownames(otu_df) <- otu_df$taxa
#  otu_df <- otu_df[, !(names(otu_df) %in% c("taxa"))]
#  otu_df["Unclassified", ] <-  rowSums(seqtab) - colSums(otu_df)
#  otu_df <- (t(apply(otu_df, 1, function(x) round(x/colSums(otu_df), digits=8)*100)))

#  return(otu_df[order(rowMeans(otu_df), decreasing = TRUE), ])
#}






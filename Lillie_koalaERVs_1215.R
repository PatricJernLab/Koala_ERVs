## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

#### Characterising polymoprhic ERVs in KOALA ####

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


#### SETUP | dir - packaages ####

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 

## set working directory

setwd("/Users/orjca337/Documents/PostDoc_UUJern/Koala_retrovirus_diversity/") 

require(dplyr)
require(data.table)
require(plyr)
require(StructuralVariantAnnotation)
require(VariantAnnotation)
require(GenomicRanges) 
require(rtracklayer)
require(regioneR)
require(stringr)
require(RColorBrewer)
require(intansv)


#### KOALA REFERENCE | fai ####

# read in koala reference fai file:
ref.fai <- fread("koala_reference/GCF_002099425.1_phaCin_unsw_v4.1_genomic.fna.fai")
ref.fai <- ref.fai[,1:2]
names(ref.fai) <- c("contig", "length")
ref.fai # 1907 contigs (1906 autosomal + mt)

## read in reference file for contigs accessions from genbank and refseq 
contig.ref <- fread("koala_reference_contig_reference.csv")
contig.ref

## read in reference file for ervs found by retrotector and their locations on genbank and refseq scaffs
erv.reference <- fread("erv_library/erv_reference.csv")
head(erv.reference)
# all the retrotector hits with score over 300
erv.scaff.map.range <- toGRanges(erv.reference[,c("hitrefseq", "scaff.start", "scaff.end", "erv.id")])
erv.scaff.map.range # 991
length(unique(erv.scaff.map.range@seqnames)) # 279 scaffolds with retrotectored ervs



#### REPEATMASKER | output from Koala reference ####
## read in repeatmasker output for the koala reference genome
ref.repeatmasker <- fread("GCF_002099425.1_phaCin_unsw_v4.1_rm_eds.out")
ref.repeatmasker <- ref.repeatmasker[3:nrow(ref.repeatmasker),]
head(ref.repeatmasker)
nrow(ref.repeatmasker) # 7814628
names(ref.repeatmasker) <- c("chr", "start", "end", "rep.class", "rep.fam")
ref.repeatmasker.range <- toGRanges(ref.repeatmasker)
ref.repeatmasker.range



#### RetroTector ERV | clade designations from phylo #### 
## read in file with details of the retrotector ervs used in the phylo and which clades they cluster into 
erv.clades <- fread("./erv_library/phaCin4v1_retrotected_retroaligned_updated1022.csv")
nrow(erv.clades)
# nb: only 129; this is only based upon the pCi phylo
head(erv.clades)
unique(erv.clades$clade)
unique(erv.clades$superclade)
erv.clades$erv.id <- erv.clades$pCi_ERV


#### BAM COV | samtools depth #### 
## read in file of bam sequencing coverage
koala.cov <- fread("koala_cov1208.csv", header = T)
koala.cov


#### BLAT | load in reference.erv.candidates ####

# curated list of candidate ervs in the reference 
# from retrotector and blat inputs
reference.erv.candidates.threeclades.final.df <- fread("reference.erv.candidates.threeclades.csv")
reference.erv.candidates.threeclades.final <- toGRanges(reference.erv.candidates.threeclades.final.df[,c("seqnames", "start","end" , "erv" ,     "score"  ,  "method", "clade")])
reference.erv.candidates.threeclades.final
# 834 ranges with clade





#### RETROSEQ | load vcfs ####
## load retroseq vcfs and convert each into grange object
koala.refgen.retroseq.vcfs <- paste("./pruned1123/", list.files(path = "./pruned1123/", 
                                                                pattern = "*vcf"), sep = "")
for (i in 1:length(koala.refgen.retroseq.vcfs)){
  name <- gsub("_Retroseq", "", gsub(".SM.*", "", gsub("_pruned.*", "", basename(koala.refgen.retroseq.vcfs[i]))))
  temp.vcf <- VariantAnnotation::readVcf(koala.refgen.retroseq.vcfs[i])
  temp.granges <- rowRanges(temp.vcf)
  temp.granges$erv <- lapply(temp.vcf@info$MEINFO, `[[`, 1)
  temp.granges$erv <- gsub("HERV\\-","HERV_",temp.granges$erv)
  temp.granges$GQ <- as.numeric(geno(temp.vcf)$GQ)
  temp.granges$SP <- as.numeric(geno(temp.vcf)$SP)
  temp.granges$CLIP3 <-as.numeric(geno(temp.vcf)$CLIP3)
  temp.granges$CLIP5 <- as.numeric(geno(temp.vcf)$CLIP5)
  temp.granges$FL <- as.numeric(geno(temp.vcf)$FL)
  names(temp.granges) <- NULL
  temp.granges$ALT <- NULL
  temp.granges$FILTER <- NULL
  temp.granges$paramRangeID <- NULL
  temp.granges$sample <- name
  if (sum(as.integer(!temp.granges@seqnames %like% "NW")) != 0) {
    temp.granges$genbank.accn <- temp.granges@seqnames
    temp.granges <- merge(temp.granges, contig.ref[,c("genbank.accn","refseq.accn")], by="genbank.accn")
    temp.granges$seqnames <- as.character(temp.granges$refseq.accn)
    temp.granges <- toGRanges(temp.granges[,c("seqnames","start","end", "REF","QUAL","erv","GQ","SP","CLIP3","CLIP5","FL","sample")])
  } else {
    temp.granges <- toGRanges(temp.granges)
  }
  assign(paste(name, ".retroseq_ref.range", sep = ""), temp.granges)
  rm(temp.vcf)
  rm(temp.granges)
}


#### RETROSEQ | define candidate loci ####
## compile retroseq vcfs | filter | reduce to candidate ERV loci 

# get retroseq ranges
obj.list = as.list(ls()[sapply(mget(ls(), .GlobalEnv), is.object)])
retroseq.granges.list <- obj.list[obj.list %like% ".retroseq_ref.range" ]
retroseq.granges.list

# compile
outgrange <- get(retroseq.granges.list[[1]])
for (i in 2:length(retroseq.granges.list)){
  tempgrange <- get(retroseq.granges.list[[i]])
  outgrange <- c(outgrange, tempgrange)
}
outgrange # 117721 ranges
unique(outgrange$sample)

## filters for fl=8 clip3>=2  clip5>=2 and adaptive filter based on cov
outgrange.fl8.clips2 <- outgrange[outgrange$FL == 8 & outgrange$CLIP3 >=2 & outgrange$CLIP5 >=2,]
outgrange.fl8.clips2 # 29187 ranges
unique(outgrange.fl8.clips2) # 24952 ranges
i=1
temp <- outgrange.fl8.clips2[outgrange.fl8.clips2$sample == unique(outgrange.fl8.clips2$sample)[i]]
bam.cov <- koala.cov$cov[koala.cov$sample %like% gsub("\\_.*", "", unique(outgrange.fl8.clips2$sample)[i]) ]
temp <- temp[temp$GQ >= bam.cov/2 & temp$GQ <= 3*bam.cov]
koala.fl8.clips2.adaptiveGQfilter <- temp
for (i in 2:length(unique(outgrange.fl8.clips2$sample))){
  temp <- outgrange.fl8.clips2[outgrange.fl8.clips2$sample == unique(outgrange.fl8.clips2$sample)[i]]
  bam.cov <- koala.cov$cov[koala.cov$sample %like% gsub("\\_.*", "", unique(outgrange.fl8.clips2$sample)[i]) ]
  temp <- temp[temp$GQ >= bam.cov/2 & temp$GQ <= 3*bam.cov]
  koala.fl8.clips2.adaptiveGQfilter <- c(koala.fl8.clips2.adaptiveGQfilter, temp)
}
koala.fl8.clips2.adaptiveGQfilter # 24769 ranges


## pad ERV calls by 100 bps (50bp either side) and reduce for loci: 
koala.fl8.clips2.adaptiveGQfilter.pad100red <- GenomicRanges::reduce(resize(koala.fl8.clips2.adaptiveGQfilter, 
                                                                         width = width(koala.fl8.clips2.adaptiveGQfilter)+(100), fix = "center"))
koala.fl8.clips2.adaptiveGQfilter.pad100red 
# 16388 ranges == candidate ERV loci across all clades







#### RETROSEQ | assign ERV to candidate loci ####
# granges list
retroseq.granges.list

# prep out frame:
koala.fl8.clips2.adaptiveGQfilter.pad100red$locus <- paste(koala.fl8.clips2.adaptiveGQfilter.pad100red)
koala.fl8.clips2.adaptiveGQfilter.pad100red.out <- data.frame(locus = paste(koala.fl8.clips2.adaptiveGQfilter.pad100red))

for (i in 1:length(retroseq.granges.list)){
  sample.grange <- get(retroseq.granges.list[[i]])
  name <- unique(sample.grange$sample)
  #sample.grange <- sample.grange[sample.grange$GQ >= bam.cov & sample.grange$GQ <= 3*bam.cov & sample.grange$FL == 8 & sample.grange$CLIP3 >=2 & sample.grange$CLIP5 >=2,]
  temp.overlap.loci <- subsetByOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red, sample.grange)
  temp.overlap.samp <- subsetByOverlaps(sample.grange, temp.overlap.loci)
  temp.overlap.samp$locus <- paste(temp.overlap.samp)
  temp.overlap.samp <- temp.overlap.samp[!duplicated(temp.overlap.samp$locus),]
  # removing duplicated regions and choosing those with the highest filter values
  temp.overlap.samp.mcols <- as.data.frame(mcols(temp.overlap.samp)) %>%
    dplyr::arrange(locus, desc(FL)) %>% 
    group_by(locus) %>%
    top_n(2, abs(FL)) %>%
    top_n(1, abs(GQ)) %>%
    as.data.frame()
  temp.overlap.samp.mcols <- temp.overlap.samp.mcols[!duplicated(temp.overlap.samp.mcols$locus),]
  temp.overlap.samp.mcols$locus.erv <- paste(temp.overlap.samp.mcols$locus, temp.overlap.samp.mcols$erv, sep="_")
  temp.overlap.samp$locus.erv <- paste(temp.overlap.samp$locus, temp.overlap.samp$erv, sep="_")
  temp.overlap.samp <- temp.overlap.samp[temp.overlap.samp$locus.erv %in% temp.overlap.samp.mcols$locus.erv,]
  
  temp.overlap.loci <- subsetByOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red, temp.overlap.samp)
  temp.hits<- findOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red, temp.overlap.samp)
  erv.col <- unique(CharacterList(split(temp.overlap.samp$erv[subjectHits(temp.hits)],
                                        queryHits(temp.hits))))
  site.col <- unique(CharacterList(split(temp.overlap.samp@ranges@start[subjectHits(temp.hits)],
                                         queryHits(temp.hits))))
  mcols(temp.overlap.loci) <- DataFrame(mcols(temp.overlap.loci), erv.col, site.col)
  names(mcols(temp.overlap.loci)) <- c("locus" , paste(name,".erv", sep = ""), paste(name,".site", sep = ""))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out <- merge(koala.fl8.clips2.adaptiveGQfilter.pad100red.out, 
                                          mcols(temp.overlap.loci), by="locus", all.x=TRUE)
}


# pull out matrix of ervs
koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout <- koala.fl8.clips2.adaptiveGQfilter.pad100red.out[, grep("erv", names(koala.fl8.clips2.adaptiveGQfilter.pad100red.out))]
rownames(koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout) <- koala.fl8.clips2.adaptiveGQfilter.pad100red.out$locus
koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout.mat <-  as.matrix(koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout)

# pull out matrix of sites
koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout <- koala.fl8.clips2.adaptiveGQfilter.pad100red.out[, grep("site", names(koala.fl8.clips2.adaptiveGQfilter.pad100red.out))]
rownames(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout) <- koala.fl8.clips2.adaptiveGQfilter.pad100red.out$locus
koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout.mat <-  as.matrix(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout)


## erv assignment to loci:
koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp <- data.frame(locus = rownames(koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout), 
                                                                           erv.assigned = NA, erv.assigned.freq = NA, 
                                                                           actual.stop = NA, actual.start= NA)
for (i in 1:nrow(koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp)){
  erv.list <- unlist(base::strsplit(as.character(unlist(t(koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout.mat[i,1:11]))),"\\-"))
  erv.table <- as.data.frame(table(unlist(erv.list), exclude = c("hybrid", "unknown", NA)))
  clade.table <- merge(data.frame(erv.id = erv.table$Var1, Freq=erv.table$Freq), erv.clades[,c("erv.id", "clade")], by="erv.id", all.x=T)
  clade.table$clade[is.na(clade.table$clade)] <- as.character(clade.table$erv.id[is.na(clade.table$clade)])
  clade.table <- aggregate(clade.table$Freq ~ clade.table$clade, clade.table, sum)
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned[i] <- paste(erv.table$Var1[erv.table$Freq == max(erv.table$Freq)], collapse="|")
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned.freq[i] <- max(erv.table$Freq) / sum(erv.table$Freq)
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$actual.stop[i] <- max(na.omit(as.numeric(unlist(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout.mat[i,1:11]))))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$actual.start[i] <- min(na.omit(as.numeric(unlist(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout.mat[i,1:11]))))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.ervs[i] <- paste(unlist(erv.table$Var1), collapse="|")
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.erv.counts[i] <- paste(unlist(erv.table$Freq), collapse="|")
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.clades[i] <- paste(unlist(clade.table$`clade.table$clade`), collapse="|")
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.clades.count[i] <- paste(unlist(clade.table$`clade.table$Freq`), collapse="|")
  ervs <- unlist(base::strsplit(koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned[i], split = "\\|"))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned.clades[i] <- ifelse(length(unique(erv.clades$clade[erv.clades$erv.id %in% ervs])) == 0, 
                                                                                                koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned[i], 
                                                                                                paste(unique(erv.clades$clade[erv.clades$erv.id %in% ervs]), collapse="|"))
  unique.ervs <- unlist(base::strsplit(koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.ervs[i], split = "\\|"))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned.ed[i] <- ifelse(length(unique(erv.clades$clade[erv.clades$erv.id %in% unique.ervs])) == 1 & length(ervs) != 1, 
                                                                                            ervs[1], koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned[i])
}





#### RETROSEQ | get counts for candidate loci #### 

koala.fl8.clips2.adaptiveGQfilter.pad100red

runs <- data.frame(FL = c(3, 3), 
                   GQ = c(1, 3))
runs

for (n in 1:nrow(runs)) {
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts <- koala.fl8.clips2.adaptiveGQfilter.pad100red
  for (i in 1:length(retroseq.granges.list)){
    g1 <- get(retroseq.granges.list[[i]])
    g1 <- g1[g1$GQ >= runs$GQ[n] & g1$FL >= runs$FL[n],]
    sample <- unique(g1$sample)
    koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts$hits <- countOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts, g1, type = "any")
    koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts$hits[koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts$hits >= 1] <- 1
    names(mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts))[names(mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts)) == "hits"] <- sample
  }
  assign(paste("koala.granges.fl8.clips2.pad100red.counts.gq", runs$GQ[n], "fl" , runs$FL[n], sep = ""), 
         koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts)
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.df <- as.data.frame(mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus <- data.frame(locus = paste(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts), 
                                                                                counts = rowSums(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.df[,2:12]))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$chr <- gsub("\\:.*", "", koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$locus)
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$start <- ranges(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts)@start
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$length <- ranges(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts)@width
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$end <- koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$start + koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$length -1
  assign(paste("koala.granges.fl8.clips2.pad100red.counts.gq", runs$GQ[n], "fl" , runs$FL[n],  ".perlocus", sep = ""), 
         koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus)
}



#### RETROSEQ | candidate loci overlaps retrotector ERVs ####

koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE <- resize(koala.fl8.clips2.adaptiveGQfilter.pad100red, 
                                                               width = width(koala.fl8.clips2.adaptiveGQfilter.pad100red)-(100), fix = "center")

# check overlaps to retrotector: erv.scaff.map.range
koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.referv.overlap <- subsetByOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE, erv.scaff.map.range)
temp.hits<- findOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE, erv.scaff.map.range)
ref.erv <- unique(CharacterList(split(erv.scaff.map.range$erv.id[subjectHits(temp.hits)],
                                      queryHits(temp.hits))))
mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.referv.overlap) <- DataFrame(mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.referv.overlap), ref.erv)
# remove these from the downstream dataset: koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.referv.overlap$locus


#### RETROSEQ | candidate loci overlaps repeatmasker ####
# check overlaps to repeatmasker: ref.repeatmasker.range
ref.repeatmasker.range.simlow <- ref.repeatmasker.range[ref.repeatmasker.range$rep.fam %in% c("Low_complexity", "Simple_repeat")]
ref.repeatmasker.range.simlow <- resize(ref.repeatmasker.range.simlow, 
                                        width = width(ref.repeatmasker.range.simlow)+(50), fix = "center")

## simple repeat overlaps: 
koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.simlow.overlaps <- subsetByOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE, ref.repeatmasker.range.simlow)
temp.hits<- findOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE, ref.repeatmasker.range.simlow)
rep.class <- unique(CharacterList(split(ref.repeatmasker.range.simlow$rep.class[subjectHits(temp.hits)],
                                        queryHits(temp.hits))))
rep.fam <- unique(CharacterList(split(ref.repeatmasker.range.simlow$rep.fam[subjectHits(temp.hits)],
                                      queryHits(temp.hits))))
mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.simlow.overlaps) <- DataFrame(mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.simlow.overlaps), rep.class, rep.fam)
koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.simlow.overlaps$rep.overlap <- 1
# remove these regions from downstream analysis: koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.simlow.overlaps$locus

# merge datasets to get info on repeats: 
koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.wrepeats <- merge(koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp, 
                                                                            mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.simlow.overlaps), 
                                                                            by="locus", all.x=TRUE)
koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.wrepeats$rep.overlap[is.na(koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.wrepeats$rep.overlap)] <- 0

# table of ERVs that overlap repeats:
koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.wrepeats.table <- as.data.frame(table(koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.wrepeats$erv.assigned.ed[koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.wrepeats$rep.overlap == 1]))

# table of ERVs that do not overlap repeats: 
koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.notrepeats.table <- as.data.frame(table(koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.wrepeats$erv.assigned.ed[koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.wrepeats$rep.overlap == 0]))


#### RETROSEQ | filtering out loci that overlap retroTector ERVs and repeatmasker simple repeat / low complexity regions ####
# remove these from the dataset
# koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.referv.overlap$locus
# length(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.referv.overlap$locus) # 149

# and these overlaps to simple repeats and low complexity regions
# koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.simlow.overlaps$locus
# length(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.simlow.overlaps$locus) # 10284

locitoremove.list <- unique(c(koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.referv.overlap$locus, 
                              koala.fl8.clips2.adaptiveGQfilter.pad100red.REALSIZE.simlow.overlaps$locus))


koala.fl8.clips2.adaptiveGQfilter.pad100red.fil <- koala.fl8.clips2.adaptiveGQfilter.pad100red[!koala.fl8.clips2.adaptiveGQfilter.pad100red$locus %in% locitoremove.list,]
length(koala.fl8.clips2.adaptiveGQfilter.pad100red.fil) 
# 5965 loci left after removing all retrotector and simple repeat/low complex overlaps


## geno matrix: 
koala.granges.fl8.clips2.pad100red.counts.gq1fl3$locus <- paste(koala.granges.fl8.clips2.pad100red.counts.gq1fl3)
koala.granges.fl8.clips2.pad100red.counts.gq1fl3.mat <- as.matrix(as.data.frame(mcols(koala.granges.fl8.clips2.pad100red.counts.gq1fl3)[2:12]))
rownames(koala.granges.fl8.clips2.pad100red.counts.gq1fl3.mat) <- koala.granges.fl8.clips2.pad100red.counts.gq1fl3$locus

koala.granges.fl8.clips2.pad100red.counts.gq1fl3.perlocus.fil <- koala.granges.fl8.clips2.pad100red.counts.gq1fl3.perlocus[koala.granges.fl8.clips2.pad100red.counts.gq1fl3.perlocus$locus %in% koala.fl8.clips2.adaptiveGQfilter.pad100red.fil$locus,]
koala.granges.fl8.clips2.pad100red.counts.gq1fl3.fil <- koala.granges.fl8.clips2.pad100red.counts.gq1fl3[koala.granges.fl8.clips2.pad100red.counts.gq1fl3$locus %in% koala.fl8.clips2.adaptiveGQfilter.pad100red.fil$locus,]


## erv assignments: 
koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.fil <- koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp[!koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$locus %in% locitoremove.list,]


#### RETROSEQ | filtered dataset | focus on three clades ####
three.erv.clades <- erv.clades[erv.clades$clade %in% c("Betalike_crown", "SMRVlike", "KoRVlike")]
koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.threeclades.fil <- koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.fil[koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.fil$erv.assigned.ed %in% three.erv.clades$pCi_ERV,]
koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat <- koala.granges.fl8.clips2.pad100red.counts.gq1fl3.mat[rownames(koala.granges.fl8.clips2.pad100red.counts.gq1fl3.mat) %in% 
                                                                                                                                                             gsub("RETROSEQ_", "", koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.threeclades.fil$locus),]

#### DELETIONS | ####

## using intansv to import data to genomic ranges
## using length cutoffs 100 and 30000 bp 
## Have noticed that Birke has many genome blocks without coverage.. genomic rearrangements in tissue? 


#### DELETIONS | load lumpy results ####
lumpy.dir <- "./lumpy"
lumpy.vcfs <- list.files(path = lumpy.dir, 
                         pattern = "*.vcf", full.names = T)

for (i in 1:length(lumpy.vcfs)){
  name <- gsub(".SM.lumpy.vcf", ".lumpy", basename(lumpy.vcfs[i]))
  loadin <- VariantAnnotation::readVcf(lumpy.vcfs[i])
  assign(name, loadin)
  loadin.intansv <- readLumpy(lumpy.vcfs[i], regSizeLowerCutoff=100, regSizeUpperCutoff=30000,method="Lumpy")
  loadin.intansv$inv <- NULL
  loadin.intansv$dup <- NULL
  assign(paste(name, ".intansv", sep = ""), loadin.intansv)
}

#### DELETIONS | load delly results ####

delly.intansv <- readDelly(file="./delly/0616_fullfullgeno_merged.vcf", 
                           regSizeLowerCutoff=100, regSizeUpperCutoff=30000,
                           readsSupport=3, method="Delly")

# also read in delly rowgranges - use chr and start to merge to intansv so we can retain DEL info
delly.vcf <- VariantAnnotation::readVcf("./delly/0616_fullfullgeno_merged.vcf") ## multiple
delly.rowgranges <- rowRanges(delly.vcf) # 115007
delly.GT <- geno(delly.vcf)$GT

DEL.ref <- data.frame(chr=delly.rowgranges@seqnames, start=delly.rowgranges@ranges@start, DEL=names(delly.rowgranges))
delly.intansv.df <- as.data.frame(delly.intansv$del)
names(delly.intansv.df) <- c("chr", "start", "end", "size", "SU")
delly.intansv.df$SU <- as.numeric(gsub("SU=", "", delly.intansv.df$SU))
delly.intansv.df <- merge(delly.intansv.df, DEL.ref, by=c("chr", "start"))
delly.intansv.df <- delly.intansv.df[order(delly.intansv.df$chr),]
delly.intansv.df[duplicated(delly.intansv.df$chrstart),]
delly.intansv.ranges <- toGRanges(delly.intansv.df[,c("chr", "start", "end")])
mcols(delly.intansv.ranges) <- delly.intansv.df[,c("SU","DEL")]
delly.intansv.ranges$locus <- paste(delly.intansv.ranges)
delly.intansv.GT <- delly.GT[rownames(delly.GT) %in% delly.intansv.ranges$DEL,]


#### DELETIONS | find overlaps to focal three clades #### 
# have: delly and lumpys 
# and have: reference.erv.candidates.threeclades.final

# can togrange instanv
list.list = as.list(ls()[sapply(mget(ls(), .GlobalEnv), is.list)])
lumpy.intansv.list <- unlist(list.list[list.list %like% "lumpy" & 
                                         list.list %like% "intansv" &
                                         !list.list %like% "list"])
lumpy.intansv.list

lumpy.summary <- data.frame(samp = gsub("\\..*", "", lumpy.intansv.list), 
                            lumpy.dels = NA, 
                            lumpy.erv.dels=NA)
for (i in 1:length(lumpy.intansv.list)){
  intansv.temp <- get(lumpy.intansv.list[i])
  name <- gsub("\\..*", "", lumpy.intansv.list[i])
  grange.temp <- toGRanges(intansv.temp$del)
  lumpy.summary$lumpy.dels[i] <- length(grange.temp)
  grange.temp.erv <- subsetByOverlaps(grange.temp, reference.erv.candidates.threeclades.final)
  if (length(grange.temp.erv) != 0){
    mcols(grange.temp.erv) <- NULL
    grange.temp.erv$method <- "lumpy"
    grange.temp.erv$sample <- name
  }
  assign(paste(name, ".lumpy.erv.ranges", sep=""), grange.temp.erv)
  lumpy.summary$lumpy.erv.dels[i] <- length(grange.temp.erv)
}
lumpy.summary

obj.list <- as.list(ls()[sapply(mget(ls(), .GlobalEnv), is.object)])
lumpy.erv.list <- obj.list[obj.list %like% ".lumpy.erv.ranges" & 
                             !obj.list %like% "all.lumpy.erv.ranges"]
lumpy.erv.list


all.lumpy.erv.ranges <- NULL
for (i in 1:length(lumpy.erv.list)){
  all.lumpy.erv.ranges <- c(get(unlist(lumpy.erv.list)[i]), all.lumpy.erv.ranges)
}
all.lumpy.erv.ranges # 215 ranges
all.lumpy.erv.ranges$size <- NULL
all.lumpy.erv.ranges$info <- NULL


all.delly.ranges <- delly.intansv.ranges # 7902 ranges
all.delly.erv.ranges <- subsetByOverlaps(all.delly.ranges, reference.erv.candidates.threeclades.final) # 46 ranges 
names(all.delly.erv.ranges) <- paste(all.delly.erv.ranges$DEL)
mcols(all.delly.erv.ranges) <- NULL
all.delly.erv.ranges$sample <- "dellyall"
all.delly.erv.ranges$method <- "delly"



#### DELETIONS | join delly and lumpy calls ####
joined.dels.threeclades.erv.ranges.winfo <- c(all.lumpy.erv.ranges,
                                              all.delly.erv.ranges)
joined.dels.threeclades.erv.ranges.winfo.red <- GenomicRanges::reduce(joined.dels.threeclades.erv.ranges.winfo)
joined.dels.threeclades.erv.ranges.winfo.red # 49 ranges 

## loop and check that the min and max lengths within reduced ranges are similar 
## i.e. both softwares are calling basically the same erv
for (i in 1:length(joined.dels.threeclades.erv.ranges.winfo.red)){
  temp <- subsetByOverlaps(joined.dels.threeclades.erv.ranges.winfo, joined.dels.threeclades.erv.ranges.winfo.red[i])
  joined.dels.threeclades.erv.ranges.winfo.red$counts[i] <- length(temp)
  joined.dels.threeclades.erv.ranges.winfo.red$min.len[i] <- min(temp@ranges@width)
  joined.dels.threeclades.erv.ranges.winfo.red$max.len[i] <- max(temp@ranges@width)
  joined.dels.threeclades.erv.ranges.winfo.red$samples[i] <- paste(temp$sample, collapse = "|")
  joined.dels.threeclades.erv.ranges.winfo.red$methods[i] <- paste(unique(temp$method), collapse = "|")
  temp.erv <- subsetByOverlaps(reference.erv.candidates.threeclades.final, joined.dels.threeclades.erv.ranges.winfo.red[i])
  joined.dels.threeclades.erv.ranges.winfo.red$erv[i] <- paste(unique(temp.erv$erv), collapse = "|")
  joined.dels.threeclades.erv.ranges.winfo.red$clade[i] <- paste(unique(temp.erv$clade), collapse = "|")
}
joined.dels.threeclades.erv.ranges.winfo.red

plot(joined.dels.threeclades.erv.ranges.winfo.red$min.len, joined.dels.threeclades.erv.ranges.winfo.red$max.len)
# good match so continue with joined.dels.threeclades.erv.ranges.winfo.red


#### DELETIONS | threeclades loci defined | get counts  #### 
# joined.dels.threeclades.erv.ranges.winfo.red  # 49 ranges



#### DELETIONS | threeclades | Delly counts #### 
# rename bilbo in delly.intansv.GT - its still in ERR name
colnames(delly.intansv.GT)[which(colnames(delly.intansv.GT) == "ERR3485163")] <- "Bilbo"
# fix it for ervs:
delly.erv.GT.reference <- delly.intansv.GT[rownames(delly.intansv.GT) %in% names(all.delly.erv.ranges),]
delly.erv.GT <- delly.intansv.GT[rownames(delly.intansv.GT) %in% names(all.delly.erv.ranges),]

# set sample names order from that in retroseq
sample.names <- colnames(koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat)
sample.names 
## use these sample names throughout

# edit genotypes to reflect ERV presence as either homo or het (1) or not present (0)
delly.erv.GT[delly.erv.GT == "0/0"] <- 1 # no deletion? Therefore has ref erv
delly.erv.GT[delly.erv.GT == "0/1"] <- 1 # heterozygote - has one copy of ref erv
delly.erv.GT[delly.erv.GT == "1/1"] <- 0 # homozygote for deletion - therefore does not have the reference erv
delly.erv.GT[delly.erv.GT == "./."] <- NA
delly.erv.GT.mat <- matrix(as.numeric(delly.erv.GT), nrow = 46, ncol = 11, byrow = F)
rownames(delly.erv.GT.mat) <- rownames(delly.erv.GT)
colnames(delly.erv.GT.mat) <- colnames(delly.erv.GT)
# check that genos are same:
head(delly.erv.GT)
head(delly.erv.GT.mat)

# join delly.erv.GT.mat to the grange
delly.erv.loci.counts <- all.delly.erv.ranges
mcols(delly.erv.loci.counts) <- delly.erv.GT.mat



#### DELETIONS | threeclades | Lumpy counts #### 
lumpy.erv.loci <- joined.dels.threeclades.erv.ranges.winfo.red[joined.dels.threeclades.erv.ranges.winfo.red$methods %like% "lumpy",] # 41 ranges
lumpy.erv.loci.counts <- lumpy.erv.loci
mcols(lumpy.erv.loci.counts) <- NULL
# getting counts - evals <= 1 
for (i in 1:length(sample.names)){
  g1 <- all.lumpy.erv.ranges[all.lumpy.erv.ranges$sample == sample.names[i]]
  lumpy.erv.loci.counts$hits <- countOverlaps(lumpy.erv.loci.counts, g1, type = "any")
  names(mcols(lumpy.erv.loci.counts))[names(mcols(lumpy.erv.loci.counts)) == "hits"] <- sample.names[i]
}
lumpy.erv.loci.counts



#### DELETIONS | threeclades | collect counts #### 

delly.erv.loci.counts
names(delly.erv.loci.counts) <- paste(names(delly.erv.loci.counts), paste(delly.erv.loci.counts), sep="_")
lumpy.erv.loci.counts
names(lumpy.erv.loci.counts) <- paste("lumpy", paste(lumpy.erv.loci.counts), sep="_")
# invert lumpy calls to reflect erv presence (no deletion) or erv absence (deletion)
lumpy.erv.loci.counts.inv <- lumpy.erv.loci.counts 
mcols(lumpy.erv.loci.counts.inv) <- sapply(mcols(lumpy.erv.loci.counts)[,1:11], function(x) 1-x)

## going to take delly genos as best caller of genotypes (eg hets), then let lumpy call the remainder - only 3 loci.
joined.erv.loci.counts <- c(delly.erv.loci.counts, lumpy.erv.loci.counts.inv[!lumpy.erv.loci.counts.inv %in% subsetByOverlaps(lumpy.erv.loci.counts.inv, delly.erv.loci.counts) ])
joined.erv.loci.counts
# 49 ranges


#### DELETION | manual corrections #### 

# manual correction after some eyeballing IGV
joined.erv.loci.counts$PacificChocolate[names(joined.erv.loci.counts) == "lumpy_NW_018344007.1:10556959-10565724"] <- 0
joined.erv.loci.counts
# remove this locus cos it aint there
joined.erv.loci.counts[names(joined.erv.loci.counts) == "lumpy_NW_018343983.1:15670484-15684950"] <- NULL
# fix bertha DEL00061256_NW_018344052.1:2620770-2644053
joined.erv.loci.counts$Bertha[names(joined.erv.loci.counts) == "DEL00061256_NW_018344052.1:2620770-2644053"] <- 1
as.data.frame(mcols(joined.erv.loci.counts))

joined.erv.loci.counts # 48 ranges




#### DELETIONS | joined.erv.loci.counts #### 
joined.dels.threeclades.erv.ranges.winfo.red.fil <- subsetByOverlaps(joined.dels.threeclades.erv.ranges.winfo.red, joined.erv.loci.counts)
joined.dels.threeclades.erv.ranges.winfo.red.fil.counts <- joined.dels.threeclades.erv.ranges.winfo.red.fil
mcols(joined.dels.threeclades.erv.ranges.winfo.red.fil.counts) <- mcols(joined.erv.loci.counts) # TEMP!! run loop to get right geno 
joined.dels.threeclades.erv.ranges.winfo.red.fil.counts
joined.dels.threeclades.erv.ranges.summary <- data.frame(locus = paste(joined.dels.threeclades.erv.ranges.winfo.red.fil.counts),
                                                         length.from.counts = NA, 
                                                         length.from.here = joined.dels.threeclades.erv.ranges.winfo.red.fil.counts@ranges@width)
for (i in 1:length(joined.dels.threeclades.erv.ranges.winfo.red.fil.counts)){
  temp <- subsetByOverlaps(joined.erv.loci.counts, joined.dels.threeclades.erv.ranges.winfo.red.fil.counts[i])
  mcols(joined.dels.threeclades.erv.ranges.winfo.red.fil.counts)[i,] <- mcols(temp)
  joined.dels.threeclades.erv.ranges.summary$length.from.counts[i] <- temp@ranges@width
}
joined.dels.threeclades.erv.ranges.summary


#### JOINED | overlaps between deletions and Retroseq ####
# deletions three clades:
joined.dels.threeclades.erv.ranges.winfo.red.fil.counts <- subsetByOverlaps(joined.dels.threeclades.erv.ranges.winfo.red.fil.counts, joined.dels.threeclades.erv.ranges.winfo.red.fil)
joined.dels.threeclades.erv.ranges.winfo.red.fil

# insertions three clades:
koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.threeclades.fil$locus <- paste("RETROSEQ", koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.threeclades.fil$locus, sep="_")
koala.granges.fl8.clips2.adaptiveGQfilter.pad100red.counts.gq1fl3.threeclades.mat <- koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat[rownames(koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat) %in% 
                                                                                                                                                        gsub("RETROSEQ_", "", koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.threeclades.fil$locus),]
# join! 
threeclades.ervdetails.retroseq <- data.frame(locus=koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.threeclades.fil$locus, 
                                              erv.assigned=koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.threeclades.fil$erv.assigned.ed, 
                                              clade=koala.fl8.clips2.adaptiveGQfilter.pad100red.ervs.threeclades.fil$erv.assigned.clades)
threeclades.ervdetails.deletions <- data.frame(locus=paste("DELETION", paste(joined.dels.threeclades.erv.ranges.winfo.red.fil), sep="_"), 
                                               erv.assigned=joined.dels.threeclades.erv.ranges.winfo.red.fil$erv, 
                                               clade=joined.dels.threeclades.erv.ranges.winfo.red.fil$clade)
threeclades.ervdetails <- rbind(threeclades.ervdetails.retroseq, threeclades.ervdetails.deletions)

## need matrix for all individuals and loci
joined.dels.threeclades.erv.ranges.winfo.red.fil.counts.df <- as.data.frame(mcols(joined.dels.threeclades.erv.ranges.winfo.red.fil.counts))
rownames(joined.dels.threeclades.erv.ranges.winfo.red.fil.counts.df) <- paste("DELETION", paste(joined.dels.threeclades.erv.ranges.winfo.red.fil), sep="_")
koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat.tojoin <- koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat
colnames(koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat.tojoin) <- gsub("\\_.*", "", colnames(koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat.tojoin))
rownames(koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat.tojoin) <- paste("RETROSEQ", rownames(koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat.tojoin), sep="_")

joined.dels.threeclades.erv.ranges.winfo.red.fil.counts.df.tojoin <- joined.dels.threeclades.erv.ranges.winfo.red.fil.counts.df[,c(colnames(koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat.tojoin))]

threeclades.ervcounts <- rbind(koala.granges.fl8.clips2.pad100red.counts.gq1fl3.threeclades.mat.tojoin, 
                               joined.dels.threeclades.erv.ranges.winfo.red.fil.counts.df.tojoin )

#### softclip ####
softclip.files <- list.files(path = "./pruned1123", pattern = "*.summary", full.names = T)
softclip.files

softclip.summary <- data.frame(samp=gsub("\\_.*", "", gsub("\\..*", "", basename(softclip.files))), 
                               locus.hits=NA)

softclip.allsamps <- NULL
for (i in 1:length(softclip.files)){
  name <-  gsub("\\_.*", "", gsub("\\..*", "", basename(softclip.files[i])))
  loadin <- fread(softclip.files[i], sep='\t')
  loadin$V5 <- NULL
  names(loadin) <- c("locus", "read", "erv.id", "match")
  loadin$samp <- name
  softclip.allsamps <- rbind(softclip.allsamps, loadin)
  softclip.summary$locus.hits[softclip.summary$samp == name] <- length(unique(loadin$locus))
}


## make a grange of the locations
uniq.softclip.loci <- data.frame(loci = unique(softclip.allsamps$locus))
uniq.softclip.loci$chr <- gsub("\\:.*", "", gsub("RETROSEQ_", "", uniq.softclip.loci$loci))
uniq.softclip.loci$start <- gsub("\\-.*",   "", gsub(".*\\:", "", gsub("RETROSEQ_", "", uniq.softclip.loci$loci)))
uniq.softclip.loci$end <- gsub(".*\\-",  "", uniq.softclip.loci$loci)
uniq.softclip.loci.grange <- toGRanges(uniq.softclip.loci[,c("chr", "start", "end", "loci")])

# find overlaps to retroseq loci
retroseq.loci <- data.frame(loci = threeclades.ervdetails$locus[threeclades.ervdetails$locus %like% "RETRO"])
retroseq.loci$chr <- gsub("\\:.*", "", gsub("RETROSEQ_", "", retroseq.loci$loci))
retroseq.loci$start <- gsub("\\-.*",   "", gsub(".*\\:", "", gsub("RETROSEQ_", "", retroseq.loci$loci)))
retroseq.loci$end <- gsub(".*\\-",  "", retroseq.loci$loci)
retroseq.loci.grange <- toGRanges(retroseq.loci[,c("chr", "start", "end", "loci")])

uniq.softclip.loci.grange.retain <- subsetByOverlaps(uniq.softclip.loci.grange, retroseq.loci.grange) # 263 ranges
hits <- findOverlaps(uniq.softclip.loci.grange, retroseq.loci.grange)
retroseq.locus <- CharacterList(split(retroseq.loci.grange$loci[subjectHits(hits)],
                               queryHits(hits)))
mcols(uniq.softclip.loci.grange.retain) <- DataFrame(mcols(uniq.softclip.loci.grange.retain), retroseq.locus)
uniq.softclip.loci.grange.retain$match <- ifelse(uniq.softclip.loci.grange.retain$loci == uniq.softclip.loci.grange.retain$retroseq.locus, 1,0)
uniq.softclip.loci.grange.retain[uniq.softclip.loci.grange.retain$match == 0]
# fix these that fo not match

softclip.allsamps.fil <- softclip.allsamps[softclip.allsamps$locus %in% uniq.softclip.loci.grange.retain$loci,]
# fix RETROSEQ_NW_018343976.1:12693058-12693160 RETROSEQ_NW_018343976.1:12693058-12693164
softclip.allsamps.fil$locus[softclip.allsamps.fil$locus == "RETROSEQ_NW_018343976.1:12693058-12693160"] <- "RETROSEQ_NW_018343976.1:12693058-12693164"
# fix RETROSEQ_NW_018343997.1:11090634-11090792 RETROSEQ_NW_018343997.1:11090687-11090792 
softclip.allsamps.fil$locus[softclip.allsamps.fil$locus == "RETROSEQ_NW_018343997.1:11090634-11090792"] <- "RETROSEQ_NW_018343997.1:11090687-11090792"

# get counts from softclip 
softclip.allsamps.fil.counts <- data.frame(locus=unique(softclip.allsamps.fil$locus))
softclip.allsamps.fil.counts.summary <- data.frame(locus=unique(softclip.allsamps.fil$locus))
sample.names <- unique(softclip.allsamps.fil$samp)
for (j in 1:length(sample.names)){
  softclip.allsamps.fil.counts[, paste(sample.names[j], "hit", sep=".")] <- NA
  softclip.allsamps.fil.counts.summary[,c(paste(sample.names[j], "ervs", sep="."))] <- NA
  softclip.allsamps.fil.counts.summary[,c(paste(sample.names[j], "ervcounts", sep="."))] <- NA
  softclip.allsamps.fil.counts.summary[,c(paste(sample.names[j], "toperv.sims", sep="."))]<- NA
  softclip.allsamps.fil.counts.summary[,c(paste(sample.names[j], "erv.clade", sep="."))]<- NA
}


for (i in 1:length(unique(softclip.allsamps.fil$locus))){
  for (j in 1:length(sample.names)){
    temp <- softclip.allsamps.fil[softclip.allsamps.fil$locus == unique(softclip.allsamps.fil$locus)[i] &
                                softclip.allsamps.fil$samp == sample.names[j] , ] %>%
      dplyr::arrange(read, desc(match)) %>% 
      group_by(read) %>%
      top_n(1, abs(match)) %>%
      dplyr::slice(1) %>%
      as.data.frame() 
    temp
    temp <- merge(temp, erv.clades[,c("erv.id", "clade")], by = "erv.id")
    # filter temp for the clade the locus was assigned to
    temp.fil.clade <- temp[temp$clade == threeclades.ervdetails$clade[threeclades.ervdetails$locus == unique(softclip.allsamps.fil$locus)[i]],]
    temp.clade <- names(sort(table(temp$clade), decreasing = TRUE)) # all clade hits
    softclip.allsamps.fil.counts[, paste(sample.names[j], "hit", sep=".")][i] <-  ifelse(nrow(temp.fil.clade)>0, 1, 0) ## only take counts for ind with softclips matching the erv clade assigned
    softclip.allsamps.fil.counts.summary[,c(paste(sample.names[j], "ervs", sep="."))][i] <- ifelse(nrow(temp)>0,
                                                                                               paste(names(sort(table(temp$erv.id), decreasing = TRUE)), collapse = "|"), NA)
    softclip.allsamps.fil.counts.summary[,c(paste(sample.names[j], "ervcounts", sep="."))][i] <- ifelse(nrow(temp)>0,
                                                                                                    paste(sort(table(temp$erv.id), decreasing = TRUE), collapse = "|"), NA)
    softclip.allsamps.fil.counts.summary[,c(paste(sample.names[j], "toperv.sims", sep="."))][i] <- ifelse(nrow(temp)>0,
                                                                                                      paste(unique(temp$erv.id[temp$match == max(temp$match)]), collapse = "|"), NA)
    softclip.allsamps.fil.counts.summary[,c(paste(sample.names[j], "erv.clade", sep="."))][i] <- ifelse(nrow(temp)>0,
                                                                                                    paste(unique(temp$clade[temp$match == max(temp$match)]), collapse = "|"), NA)
  }
  softclip.allsamps.fil.counts.summary$uniqueervs[i] <- paste(unique(na.omit(unlist(strsplit(as.character(softclip.allsamps.fil.counts.summary[,grep("\\.ervs", names(softclip.allsamps.fil.counts.summary))][i,]), "\\|")))), collapse = "|")
  softclip.allsamps.fil.counts.summary$uniqueclades[i] <- paste(unique(na.omit(unlist(strsplit(as.character(softclip.allsamps.fil.counts.summary[,grep("erv.clade", names(softclip.allsamps.fil.counts.summary))][i,]), "\\|")))), collapse = "|")
}


#### softclips rbind counts from softclip supported locus to 
rownames(softclip.allsamps.fil.counts) <- softclip.allsamps.fil.counts$locus
softclip.allsamps.fil.counts$locus <- NULL
colnames(softclip.allsamps.fil.counts) <- gsub("\\.hit", "", colnames(softclip.allsamps.fil.counts))
threeclades.ervcounts.remainder <- threeclades.ervcounts[!rownames(threeclades.ervcounts) %in% rownames(softclip.allsamps.fil.counts),]
softclip.allsamps.fil.loci <- rownames(softclip.allsamps.fil.counts)

threeclades.ervcounts.softclips <- rbind(softclip.allsamps.fil.counts, threeclades.ervcounts.remainder)


threeclades.ervdetails.softclips <- threeclades.ervdetails
threeclades.ervdetails.softclips$locus <- as.character(threeclades.ervdetails.softclips$locus)
threeclades.ervdetails.softclips$locus[threeclades.ervdetails.softclips$locus %in% softclip.allsamps.fil.loci] <- paste(threeclades.ervdetails.softclips$locus[threeclades.ervdetails.softclips$locus %in% softclip.allsamps.fil.loci], "s", sep="")



### heatmaps #### 

map.plot.details <- data.frame(clade=c("Betalike_crown", "KoRVlike", "SMRVlike"),
                               mapcolref=c(1, 5, 3))
mapcol=brewer.pal(name = "Paired", n = 6)
rename.df <- data.frame(samp.names = c("Amity","Bertha","Bilbo","Birke","Guppy","Indigo","Jaffa", "Meander","Mintie","PacificChocolate", "Utopia"), 
                        id.names=c("LC1", "LC2", "HC1","HC2","LC3","LC4","LC5","LC6","LC7","HC3","LC8"))
rename.df

colnames(threeclades.ervcounts.softclips)
threeclades.ervcounts.softclips.renamed <- threeclades.ervcounts.softclips
colnames(threeclades.ervcounts.softclips.renamed) <- rename.df$id.names
head(threeclades.ervcounts.softclips.renamed)

for (i in 1:nrow(map.plot.details)){
  genomat <- as.matrix(threeclades.ervcounts.softclips.renamed[rownames(threeclades.ervcounts.softclips.renamed) %in% gsub("s", "", threeclades.ervdetails.softclips$locus[threeclades.ervdetails.softclips$clade == as.character(map.plot.details$clade[i])]),])
  genomat <- genomat[,c("HC1","HC2","HC3","LC1","LC2","LC3","LC4","LC5","LC6","LC7","LC8")]
  genomat <- genomat[order(genomat[,c("LC8")]),]
  genomat <- genomat[order(genomat[,c("LC7")]),]
  genomat <- genomat[order(genomat[,c("LC6")]),]
  genomat <- genomat[order(genomat[,c("LC5")]),]
  genomat <- genomat[order(genomat[,c("LC4")]),]
  genomat <- genomat[order(genomat[,c("LC3")]),]
  genomat <- genomat[order(genomat[,c("LC2")]),]
  genomat <- genomat[order(genomat[,c("LC1")]),]
  genomat <- genomat[order(genomat[,c("HC3")]),]
  genomat <- genomat[order(genomat[,c("HC2")]),]
  genomat <- genomat[order(genomat[,c("HC1")]),]
  pdf(paste("1215_checks_", map.plot.details$clade[i], ".pdf", sep=""))
  heatmap(genomat , scale = "none", 
          Rowv = NA, Colv = NA, cexCol=2, 
          col = mapcol[map.plot.details$mapcolref[i]:(map.plot.details$mapcolref[i]+1)], 
          main = paste(nrow(genomat), map.plot.details$clade[i], "hits"),
          margins = c(4,15)) 
  dev.off()
}


## done

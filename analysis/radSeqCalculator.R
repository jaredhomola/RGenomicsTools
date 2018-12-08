#################################### 
####     RAD-seq calculator     ####
#################################### 

## Set parameters, currently assuming HiSeq X
genomeSize <- 598663367            ## Genome size in bp (e.g., reference genome size)
cutSites <- genomeSize/4^6         ## Number of enzyme cut sites (If PstI use 4^6, if SbfI use 4^8)
radTags <- cutSites*2              ## Number of RAD tags created from enzyme cuts
readLength <- 150*2                ## Read length (e.g., 150 for 150bp). Include *2 if paired reads
basesCovered <- radTags*readLength ## Number of bases to sequence based on read length (if 2x150 use 300)
laneReads <- 400000000             ## Reads per lane minus PhiX spike-in (i.e., 450000000 for HiSeq X, minus 50000000 for PhiX)
qualityFilter <- 0.8               ## Proportion of reads retained after quality filter
dupRemoval <- 0.75                 ## Proportion of reads retained after removing PCR duplicates
desiredDepth <- 20                 ## Desired sequencing depth

### Calculate proportion of genome to be sequenced
basesCovered/genomeSize

### Number of raw reads to sequence one individuals to desiredDepth at one locus
rawReadsNeeded.locus <- (desiredDepth / dupRemoval) / qualityFilter

### Number of raw reads to sequence per individual to acheive desiredDepth
rawReadsNeeded.ind <- rawReadsNeeded.locus * radTags

### Number of individuals that can fit on one lane to acheive desired depth at all loci
laneReads / rawReadsNeeded.ind
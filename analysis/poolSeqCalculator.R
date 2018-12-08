#########################################   
#####  Pool-seq coverage estimator  #####  
#########################################

## Set parameters, currently assuming HiSeq X
genomeSize <- 598663367            ## Genome size in bp (e.g., reference genome size)
readLength <- 150*2                ## Read length (e.g., 150 for 150bp). Include *2 if paired reads
laneReads <- 400000000             ## Reads per lane minus PhiX spike-in (i.e., 450000000 for HiSeq X, minus 50000000 for PhiX)
qualityFilter <- 0.8               ## Expected proportion of reads retained after quality filter
dupRemoval <- 0.9                  ## Expected proportion of reads retained after removing optical duplicates (common to HiSeq X)
desiredDepth <- 100                ## Desired sequencing depth


#### Assuming PCR-free library prep and HiSeq X platform
readsRetainedAfterDedup <- 0.9
desiredDepth.pooled <- 100

# Number of reads needed to generate 1x coverage across genome
rawReads.pool <- (genomeSize/readLength)

# Raw reads needed for 1x coverage following optical duplication on HiSeq X and quality filtering
rawReadsFiltered.pool <- (rawReads.pool/dupRemoval) / qualityFilter

# Raw reads needed per pool to get desiredDepth
rawReadsFiltered.pool * desiredDepth

# Number of pools per lane
laneReads / (rawReadsFiltered.pool*desiredDepth)

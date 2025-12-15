#### Libraries
library(data.table)
library(foreach)
library(poolfstat) 
library(SeqArray)

#### Data loading
meta_file <- "/scratch/cqh6wn/Class/full_sample_metadata.90Sept2025_ExpEvo.csv"
gds_file  <- "/scratch/cqh6wn/Class/dest.all.PoolSNP.001.50.29Sept2025_ExpEvo.norep.ann.gds"
out_dir   <- "/scratch/cqh6wn/Class/baypass_project_updated"

dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "results"), recursive = TRUE, showWarnings = FALSE)
### Filter metadata
samps <- fread(meta_file, fill=T)
valid_types <- c("wild", "Wild flies", "F1")
samps <- samps[Recommendation == "Pass" & fly_type %in% valid_types & set != "ExpEvo" & city != "Charlottesville"]
samps <- samps[order(jday)] # Order by date collected
candidates <- samps[, .SD[c(1, .N)], by = .(locality, year)] # group data by "locality, year", subset by selecting the first and last row
## Define spring and fall by the gap between jday
candidates[, gap := max(jday) - min(jday), by = .(locality, year)]
final_samps <- candidates[gap > 60]
final_samps[, season_defined := ifelse(jday == min(jday), "Spring", "Fall"), by = .(locality, year)]

#### Contrast
final_samps <- final_samps[order(sampleId)]
final_samps[, contrast_code := ifelse(season_defined == "Spring", 1, -1)] # If the season is Spring, assign value 1
### Save contrast
contrast_path <- file.path(out_dir, "inputs", "contrast_real.txt")
write.table(t(final_samps$contrast_code),  # Transpose
            file = contrast_path, 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#### Analyze
genofile <- seqOpen(gds_file) # gds file contains raw DNA data for all fly populations in the DEST
seqSetFilter(genofile, sample.id = final_samps$sampleId) # Focus on the final_samps flies (spring-fall pairs)
### Filter for biallelic SNPs only
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAlleles=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id"),
                     refAllele=seqGetData(genofile, "$ref"),
                     altAllele=seqGetData(genofile, "$alt"))

snp.dt <- snp.dt[nAlleles==2]
seqSetFilter(genofile, sample.id=final_samps$sampleId, variant.id=snp.dt$id)
### get allele frequencies 
ad <- seqGetData(genofile, "annotation/format/AD")$data # Number of reads for the alternative allele
dp <- seqGetData(genofile, "annotation/format/DP") # Total number of reads
#Add metadata ad
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")
#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")
ad[is.na(ad)] <- 0
dp[is.na(dp)] <- 0
## Calculate reference 
ref <- dp - ad
ref[ref < 0] <- 0 # Number of reads for the reference allele

### make  poolfstat object
snp.meta <- data.frame(Chromosome = snp.dt$chr, 
                       Position = snp.dt$pos, 
                       Ref = snp.dt$ref, 
                       Alt = snp.dt$alt)

#### Slice poolfstat Object 
num_batches <- 50

for(i in 1:num_batches) {
  index <- seq(i, nrow(snp.dt), by = num_batches)
  subpool <- new("pooldata",
                 npools = nrow(ref),              
                 nsnp = length(index),            
                 refallele.readcount = t(ref[, index]), # Transpose 
                 readcoverage = t(dp[, index]),         # Transpose 
                 poolsizes = final_samps$nFlies * 2,
                 poolnames = rownames(ref),
                 snp.info = snp.meta[index, ])
  pooldata2genobaypass(subpool, 
                       writing.dir = file.path(out_dir, "inputs"), 
                       prefix = paste0("subpool_", i))
}

seqClose(genofile)

#### Control file
job_list <- data.table(
  subpool_id = 1:num_batches,
  contrast = contrast_path, # Path of the contrast file
  output = paste0("real_run_subpool_", 1:num_batches)
)

write.table(job_list, file.path(out_dir, "baypass_control_file.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

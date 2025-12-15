#### Libraries
library(data.table)
library(SeqArray)
library(foreach)
library(doParallel)
library(ggplot2)
library(patchwork)

#### File path
gds_file    <- "/scratch/cqh6wn/Class/dest.all.PoolSNP.001.50.29Sept2025_ExpEvo.norep.ann.gds"
results_dir <- "/scratch/cqh6wn/Class/baypass_project_updated/results"
num_batches <- 50

#### SNP map
genofile <- seqOpen(gds_file)
seqResetFilter(genofile)

snp.dt <- data.table(
  chr = seqGetData(genofile, "chromosome"),
  pos = seqGetData(genofile, "position"),
  id = seqGetData(genofile, "variant.id"),
  nAlleles = seqGetData(genofile, "$num_allele")
)
snp.dt <- snp.dt[nAlleles == 2] # Flilter for biallelic SNPs

seqClose(genofile)

#### Result columns
snp.dt[, XtXst := NA_real_] # Differentiation 
snp.dt[, Contrast := NA_real_] # Seasonality

#### Aggregate results
for(i in 1:num_batches) {
  idx <- seq(i, nrow(snp.dt), by = num_batches)
  xtx_file <- file.path(results_dir, paste0("real_run_subpool_", i, "_summary_pi_xtx.out"))
  if(file.exists(xtx_file)) {
    xtx_tmp <- fread(xtx_file, header = TRUE)
    snp.dt[idx, XtXst := xtx_tmp$XtXst]
  }
  cnt_file <- file.path(results_dir, paste0("real_run_subpool_", i, "_summary_contrast.out"))
  if(file.exists(cnt_file)) {
    cnt_tmp <- fread(cnt_file, header = TRUE)
    snp.dt[idx, Contrast := abs(cnt_tmp$C2)]
  }
}

#### Filter to chromosome arms
genome_dt <- snp.dt[chr %in% c("2L", "2R", "3L", "3R") & !is.na(XtXst) & !is.na(Contrast)] # Focusing on major chromosome arms
genome_dt[, chr := factor(chr, levels = c("2L", "2R", "3L", "3R"))]
genome_dt[, q_xtx := 1 - rank(XtXst) / (.N + 1)] 
genome_dt[, q_contrast  := 1 - rank(Contrast) / (.N + 1)] # Convert raw differentiation to percentiles

#### Anova tests
anova_XtX <- aov(XtXst ~ chr, data = genome_dt) # Is XtX significantly differentiated between chromosome arms
anova_C <- aov(Contrast ~ chr, data = genome_dt)
summary(anova_XtX)
summary(anova_C)

#### Identify the most differentiated chromosome
tukey_XtX<- TukeyHSD(anova_XtX) #  Which chromosomes are different from each other
tukey_C <- TukeyHSD(anova_C)
print(tukey_XtX)
print(tukey_C)

#### Figure 2: XtX
### XtX by Chromosome
p2_a <- ggplot(genome_dt, aes(x = chr, y = XtXst, fill = chr)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = NULL, y = "XtXst") +
  theme_bw() + theme(legend.position = "none")
### In(2L)t vs. rest of the genome
genome_dt[, Region_Group := factor(fifelse(chr=="2L" & pos %between% c(2225744, 13154180), "In(2L)t",
                                           fifelse(chr=="2L", "2L", "Other Autosomes")), # Specify SNPs inside In(2L)t, SNPs on 2L but outside In(2L)t, and other genomic regions
                                   levels = c("In(2L)t", "2L", "Other Autosomes"))]
fisher.test(genome_dt$Region_Group == "In(2L)t", genome_dt$q_xtx <= 0.05) # Fisher's exact test, does In(2L)t contain top 5% differentiated SNPs at a proportion hiogher than random chance
### p2b
plot_data <- genome_dt[, .(Percent_Outliers = mean(q_xtx <= 0.05) * 100), by = Region_Group] # proportion of outliers in each region group
p2_b <- ggplot(plot_data, aes(x = Region_Group, y = Percent_Outliers, fill = Region_Group)) +
  geom_col() +
  geom_hline(yintercept = 5, linetype = "dashed", color = "black") + # the expected 5% baseline
  annotate("text", x = 3, y = 5, label = "Expected (5%)", size = 3) +
  scale_fill_manual(values = c("In(2L)t" = "red", "2L" = "steelblue", "Other Autosomes" = "gray70")) +
  labs(y = "% of SNPs that are Outliers",
       x = NULL) +
  theme_bw() +
  theme(legend.position = "none")
### Sliding window
plot_dt_c <- genome_dt[chr == "2L"] # Focusing on 2L
win_size <- 1e5; step_size <- 5e4
windows <- data.table(start = seq(min(plot_dt_c$pos), max(plot_dt_c$pos) - win_size, by = step_size))
windows[, end := start + win_size]
  
xtx_stats <- foreach(i = 1:nrow(windows), .combine = "rbind") %dopar% {
    w <- plot_dt_c[pos >= windows$start[i] & pos <= windows$end[i]] # Grabs SNPs fall inside the current window
    if(nrow(w) > 0) {
      n_out <- sum(w$q_xtx <= 0.05) # Out of the 500 SNPs in this window, how many are in the top 5% of global differentiation?
      p <- pbinom(n_out - 1, nrow(w), 0.05, lower.tail = FALSE) # What is the probability of getting greater than 0.05
      data.table(mid = (windows$start[i] + windows$end[i])/2, p = p)
    }
}
inversion_start <- 2225744
inversion_end   <- 13154180
p2_c <- ggplot(xtx_stats, aes(x = mid/1e6, y = -log10(p))) +
    geom_line(color = "black") +
    annotate("rect", xmin = inversion_start/1e6, xmax = inversion_end/1e6, ymin = 0, ymax = Inf, fill = "blue", alpha = 0.2) + # Inversion region
    labs(x = "Position (Mb)", y = "-log10(P-value)") +
    theme_bw()
### Patchwork
design <- "
  AB
  CC
"
final_figure_2 <- p2_a + p2_b + p2_c + 
  plot_layout(design = design) + 
  plot_annotation(tag_levels = 'A')

#### Figure 3: C2
### Sliding window
win_size_c2 <- 4e5; step_size <- 5e4
windows_c2 <- data.table(start = seq(min(plot_dt_c$pos), max(plot_dt_c$pos) - win_size_c2, by = step_size))
windows_c2[, end := start + win_size_c2]

cnt_stats <- foreach(i = 1:nrow(windows_c2), .combine = "rbind") %dopar% {
  w <- plot_dt_c[pos >= windows_c2$start[i] & pos <= windows_c2$end[i]]
  if(nrow(w) > 0) {
    n_out <- sum(w$q_contrast <= 0.05)
    p <- pbinom(n_out - 1, nrow(w), 0.05, lower.tail = FALSE)
    data.table(mid = (windows_c2$start[i] + windows_c2$end[i])/2, p = p)
  }
}

nunez_peaks <- data.table(peak = c(3.1, 4.7, 5.2, 6.1, 6.8, 9.6))
nunez_peaks[, `:=`(xmin=peak-0.15, xmax=peak+0.15, label=paste0(peak,"Mb"))] # 150bp flanking the peaks = the peak regions

p3_a <- ggplot() +
  geom_rect(data = nunez_peaks, aes(xmin=xmin, xmax=xmax, ymin=0, ymax=Inf), fill = "gold", alpha = 0.4) +
  geom_line(data = cnt_stats, aes(x=mid/1e6, y=-log10(p)), color = "darkred", linewidth = 0.8) +
  geom_text(data = nunez_peaks, aes(x=peak, y=13, label=label), vjust = -0.5, size = 3, fontface = "bold") +
  coord_cartesian(xlim = c(2.2, 13.2)) + # focyuing inside the inversion
  labs( x = "Position (Mb)", y = "-log10(P-value)") +
  theme_bw()

### Nunez paper enrichment analysis
nunez_peaks <- c(3.1, 4.7, 5.2, 6.1, 6.8, 9.6)
peaks_dt <- data.table(name  = paste0("Peak_", nunez_peaks),
                       start = (nunez_peaks - 0.15) * 1e6,
                       end   = (nunez_peaks + 0.15) * 1e6) # Define enrichment regions, 150bp flanking the Nunez peaks
inv_dt <- plot_dt_c[pos %between% c(2225744, 13154180)]
inv_dt[, Outlier := q_contrast <= 0.05]
### Fisher's test
overlap_stats <- peaks_dt[, {
  peak_region <- inv_dt$pos >= start & inv_dt$pos <= end # which SNPs in inv_dt fall inside the current Nunez peak
  tryCatch({
    ft <- fisher.test(inv_dt$Outlier, peak_region, alternative = "greater")  # Is this peak richer in outliers than the background?
    list(
      Outliers = sum(inv_dt$Outlier & peak_region),
      Total = sum(peak_region),
      OR = round(ft$estimate, 2),
      P_Value = ft$p.value
    )
  }, error = function(e) list(Outliers=0, Total=0, OR=NA, P_Value=NA))
}, by = name]
print(overlap_stats)
### Enrichment figure
plot_data <- data.table(
  Region = c("Peak 3.1", "Peak 4.7", "Peak 5.2", "Peak 6.1", "Peak 6.8", "Peak 9.6"),
  OR = c(0.98, 0.95, 1.09, 1.04, 1.13, 0.96),
  Sig = c("ns", "ns", "*", "ns", "***", "ns")
)
p3_b <- ggplot(plot_data, aes(x = Region, y = OR, fill = Sig != "ns")) +
  geom_col(width = 0.6, alpha = 0.9, color = "black") +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray70")) +
  geom_text(aes(label = Sig), vjust = -0.5, size = 5, fontface = "bold") +
  labs(y = "Odds Ratio", x = NULL) +
  theme_bw() +
  theme(legend.position = "none")
### Patchwork
design <- "
  AA
  BB
"
final_figure_3 <- p3_a + p3_b + 
  plot_layout(design = design) + 
  plot_annotation(tag_levels = 'A')


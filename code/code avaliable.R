
###############################################################################
# LSS and GRPS Processing Script
# 
# Prepare:
# 1. Loci information
# 2. Regulatory network
#
###############################################################################



# =============================================================================
#  Load required packages
# =============================================================================
library(data.table)
library(GenomicRanges)
library(ieugwasr)



# =============================================================================
#  Main analysis
# =============================================================================



# -----------------------------------------------------------------------------
# 1. Locus-specific stratification (LSS)
# -----------------------------------------------------------------------------

# Identify the inflection point in QQ-plots for each locus to perform LSS analysis


### Load example data
example <- read.table("LSS-and-GRPS/input/locus_example_gwas_information.txt",header = T)
length(unique(example$LeadSNP))

### Identify unique loci
cuts <- unique(example$LeadSNP)

### Initialize results data frame
infle <- as.data.frame(as.matrix(1))
infle$inflx <- 1
infle$infly <- 1
names(infle) <- c("locus","inflection_x","inflection_y")
set <- infle


### Process each locus
locus <- 1
dat <- example[example$LeadSNP == cuts[locus],]
dat <- dat[order(dat$log10P,decreasing = T ),]
obs <- dat$log10P
N <- length(obs)
theoretical <- -log(1:N/N,10)
df <- as.data.frame(cbind(theoretical,obs))

# Sort by theoretical (x) ascending
df <- df[order(df$theoretical), ]

# Fit smoothing spline
spl <- smooth.spline(df$theoretical, df$obs)

# Create uniform grid for interpolation and derivatives
x_min <- min(df$theoretical)
x_max <- max(df$theoretical)
xnew <- seq(x_min, x_max, length.out = 1000)

# Predict smoothed y, first deriv, second deriv
ynew <- predict(spl, xnew)$y
d1 <- predict(spl, xnew, deriv = 1)
d2 <- predict(spl, xnew, deriv = 2)

# Check for NA values and stop if present
if (any(is.na(d2$y))) {
  stop("Second derivative contains NA values. Check spline fit or data range.")
}

# Find indices where second deriv changes sign using d2$y
sign_changes <- which(sign(d2$y[-length(d2$y)]) != sign(d2$y[-1]))

# For each sign change, approximate x where d2$y=0, and classify
inflections <- data.frame(x = numeric(), y = numeric(), d1 = numeric(), type = character())
for (i in sign_changes) {
  # Linear interpolation for exact zero
  dx <- xnew[i+1] - xnew[i]
  zero_frac <- -d2$y[i] / (d2$y[i+1] - d2$y[i])
  x_zero <- xnew[i] + zero_frac * dx
  
  # Predict y and d1 at x_zero
  y_zero <- predict(spl, x_zero)$y
  d1_zero <- predict(spl, x_zero, deriv = 1)$y
  
  # Classify: pos to neg (local max d1, peak top), neg to pos (local min d1)
  if (d2$y[i] > 0 && d2$y[i+1] < 0) {
    inf_type <- "peak_top (convex to concave)"
  } else if (d2$y[i] < 0 && d2$y[i+1] > 0) {
    inf_type <- "valley (concave to convex)"
  } else {
    next  # Skip if no clear change
  }
  
  inflections <- rbind(inflections, data.frame(x = x_zero, y = y_zero, d1 = d1_zero, type = inf_type))
}

# Select convex bulge inflections (peak tops)
convex_bulges <- inflections[inflections$type == "peak_top (convex to concave)", ]

# Compute global max d1 from the entire first derivative
global_max_d1 <- max(d1$y)
threshold <- 0.3 * global_max_d1
filtered_bulges <- convex_bulges[convex_bulges$d1 > threshold, ]

# select the one with max x
if (nrow(filtered_bulges) > 0) {
  primary_bulge <- filtered_bulges[which.max(filtered_bulges$x), ]
  cat("Primary convex bulge inflection point (filtered and max x):\n")
  print(primary_bulge)
} else {
  cat("No convex bulge inflections found after filtering. Try adjusexampleg smoothing or threshold.\n")
}

# Plot for visualization (mark only the primary if exists)
pdf(file=paste("locus with lead snp ",cuts[locus],".pdf",sep=''),width = 4.5,height = 5)
plot(df$theoretical, df$obs, 
     pch = 19, 
     col = "gray60",          
     cex = 0.4,              
     main = cuts[locus],
     xlab = "Expected -log10 P-value", 
     ylab = "Observed -log10 P-value",
)
lines(xnew, ynew, col = "#00b4d8", lwd = 1)
if (nrow(filtered_bulges) > 0) {
  points(primary_bulge$x, primary_bulge$y, pch = 16, col = "darkred", cex = 1.5)
  text(primary_bulge$x, primary_bulge$y, labels = sprintf("(%.2f, %.2f)", round(primary_bulge$x, 2), round(primary_bulge$y, 2)), pos = 4)
}
dev.off()
graphics.off()

infle_x <- primary_bulge$x
infle_y <- primary_bulge$y
set$inflection_x <- infle_x
set$inflection_y <- infle_y
set$locus <- cuts[locus]
set$max <- max(dat$log10P)
infle <- set


for (locus in 2:length(cuts)) {
  dat <- example[example$LeadSNP == cuts[locus],]
  dat <- dat[order(dat$log10P,decreasing = T ),]
  obs <- dat$log10P
  N <- length(obs)
  theoretical <- -log(1:N/N,10)
  df <- as.data.frame(cbind(theoretical,obs))
  
  
  # Sort by theoretical (x) ascending
  df <- df[order(df$theoretical), ]
  
  # Fit smoothing spline
  spl <- smooth.spline(df$theoretical, df$obs)
  
  # Create uniform grid for interpolation and derivatives
  x_min <- min(df$theoretical)
  x_max <- max(df$theoretical)
  xnew <- seq(x_min, x_max, length.out = 1000)
  
  # Predict smoothed y, first deriv, second deriv
  ynew <- predict(spl, xnew)$y
  d1 <- predict(spl, xnew, deriv = 1)
  d2 <- predict(spl, xnew, deriv = 2)
  
  # Check for NA values and stop if present
  if (any(is.na(d2$y))) {
    stop("Second derivative contains NA values. Check spline fit or data range.")
  }
  
  # Find indices where second deriv changes sign using d2$y
  sign_changes <- which(sign(d2$y[-length(d2$y)]) != sign(d2$y[-1]))
  
  # For each sign change, approximate x where d2$y=0, and classify
  inflections <- data.frame(x = numeric(), y = numeric(), d1 = numeric(), type = character())
  for (i in sign_changes) {
    # Linear interpolation for exact zero
    dx <- xnew[i+1] - xnew[i]
    zero_frac <- -d2$y[i] / (d2$y[i+1] - d2$y[i])
    x_zero <- xnew[i] + zero_frac * dx
    
    # Predict y and d1 at x_zero
    y_zero <- predict(spl, x_zero)$y
    d1_zero <- predict(spl, x_zero, deriv = 1)$y
    
    # Classify: pos to neg (local max d1, peak top), neg to pos (local min d1)
    if (d2$y[i] > 0 && d2$y[i+1] < 0) {
      inf_type <- "peak_top (convex to concave)"
    } else if (d2$y[i] < 0 && d2$y[i+1] > 0) {
      inf_type <- "valley (concave to convex)"
    } else {
      next  # Skip if no clear change
    }
    
    inflections <- rbind(inflections, data.frame(x = x_zero, y = y_zero, d1 = d1_zero, type = inf_type))
  }
  
  # Select convex bulge inflections (peak tops)
  convex_bulges <- inflections[inflections$type == "peak_top (convex to concave)", ]
  
  # Compute global max d1 from the entire first derivative
  global_max_d1 <- max(d1$y)
  threshold <- 0.3 * global_max_d1
  filtered_bulges <- convex_bulges[convex_bulges$d1 > threshold, ]
  
  # select the one with max x
  if (nrow(filtered_bulges) > 0) {
    primary_bulge <- filtered_bulges[which.max(filtered_bulges$x), ]
    cat("Primary convex bulge inflection point (filtered and max x):\n")
    print(primary_bulge)
  } else {
    cat("No convex bulge inflections found after filtering. Try adjusexampleg smoothing or threshold.\n")
  }

  # Plot for visualization (mark only the primary if exists)
  pdf(file=paste("locus with lead snp ",cuts[locus],".pdf",sep=''),width = 4.5,height = 5)
  plot(df$theoretical, df$obs, 
       pch = 19, 
       col = "gray60",          
       cex = 0.4,               
       main = cuts[locus],
       xlab = "Expected -log10 P-value", 
       ylab = "Observed -log10 P-value",
  )
  lines(xnew, ynew, col = "#00b4d8", lwd = 1)
  if (nrow(filtered_bulges) > 0) {
    points(primary_bulge$x, primary_bulge$y, pch = 16, col = "darkred", cex = 1.5)
    text(primary_bulge$x, primary_bulge$y, labels = sprintf("(%.2f, %.2f)", round(primary_bulge$x, 2), round(primary_bulge$y, 2)), pos = 4)
  }
  dev.off()
  graphics.off()
  
  infle_x <- primary_bulge$x
  infle_y <- primary_bulge$y
  set$inflection_x <- infle_x
  set$inflection_y <- infle_y
  set$locus <- cuts[locus]
  set$max <- max(dat$log10P)
  infle <- rbind(set,infle)
}

# infle$inflection_y[infle$inflection_y < -log10(0.00000005)] <- -log10(0.00000005)
infle$pro <- infle$inflection_y/infle$max
infle$top <- 1-infle$pro




# Select variants above the threshold as candidate causal variants
r1 <- infle
r1 <- r1[,c(1,3)]
all <- example
length(unique(all$LeadSNP))
a <- unique(all$LeadSNP)

for (i in 1:length(a)) {
  i=1
  locus <- all[all$LeadSNP == a[i],]
  cutoff <- r1[r1$locus == a[i],]
  cutoff <- as.numeric(cutoff$inflection_y)
  locus_filter <- locus[locus$log10P >= cutoff,]
}

data <- locus_filter

for (i in 2:length(a)) {
  locus <- all[all$LeadSNP == a[i],]
  cutoff <- r1[r1$locus == a[i],]
  cutoff <- as.numeric(cutoff$inflection_y)
  locus_filter <- locus[locus$log10P >= cutoff,]
  data <- rbind(data,locus_filter)
}
length(unique(data$LeadSNP))
set <- data
length(unique(data$RSID))
set$Locus <- paste(set$chr_hg38_LeadSNP,set$window_left,set$window_right,sep = "_")
names(set)
set <- set[,c(1,4,25,10,21:23)]
write.csv(set,"locus_example_lss_candidate_variants.csv",row.names = F)



# -----------------------------------------------------------------------------
# 2. Functional annotation
# -----------------------------------------------------------------------------

# examined the functionality of LSS-variants

info <- read.csv("./LSS-and-GRPS/output/locus_example_lss_candidate_variants.csv",header = T)
length(unique(info$LeadSNP))
length(unique(info$RSID))
names(info)
snp <- info
length(unique(snp$RSID))

cree <- read.table("LSS-and-GRPS/input/kidney_scATACseq_peaks.txt",header = T)
snp.g=GRanges(data.frame(chrom=snp$chr_hg38,start=snp$start_hg38,end=snp$end_hg38))
cree.g=GRanges(data.frame(chrom=cree$chr_atac,start=cree$start_atac,end=cree$end_atac))
snp_cree <- suppressWarnings(findOverlaps(cree.g,snp.g))
snp_cree_overlap <- data.frame(cree[queryHits(snp_cree),],snp[subjectHits(snp_cree),])
snp_cree_overlap <-snp_cree_overlap[!duplicated(snp_cree_overlap),]
names(snp_cree_overlap)
length(unique(snp_cree_overlap$RSID))
length(unique(snp_cree_overlap$LeadSNP))
func <- snp_cree_overlap
names(func)
func <- func[,c(5:11)]
func <- func[!duplicated(func),]
write.csv(func,"./LSS-and-GRPS/output/locus_example_lss_functional_variants.csv",row.names = F)

# Identification of candidate causal genes
reg <- readRDS("./LSS-and-GRPS/input/kidney_regulation.rds")
snp <- func
library(GenomicRanges)
snp.g=GRanges(data.frame(chrom=snp$chr_hg38,start=snp$start_hg38,end=snp$end_hg38))
reg.g=GRanges(data.frame(chrom=reg$chr,start=reg$start,end=reg$end))
snp_reg <- suppressWarnings(findOverlaps(reg.g,snp.g))
snp_reg_overlap <- data.frame(reg[queryHits(snp_reg),],snp[subjectHits(snp_reg),])
snp_reg_overlap <-snp_reg_overlap[!duplicated(snp_reg_overlap),]
info_reg <- snp_reg_overlap
length(unique(info_reg$LeadSNP))
length(unique(info_reg$RSID))
length(unique(info_reg$gene))
names(info_reg)
write.csv(info_reg,"./LSS-and-GRPS/output/locus_example_lss_variants_regulated_genes.csv",row.names = F)


# -----------------------------------------------------------------------------
# 3. Gene regulatory prioritization score (GRPS)
# -----------------------------------------------------------------------------

#Take full consideration of the complexity of transcriptional regulation to prioritize genes

locus <- read.csv("./LSS-and-GRPS/output/locus_example_lss_variants_regulated_genes.csv",header = T)
length(unique(locus$gene))
length(unique(locus$RSID))
length(unique(locus$LeadSNP))
names(locus)
locus$peak <- paste(locus$chr,locus$start,locus$end,sep = '_')
gwas <- read.table("LSS-and-GRPS/input/locus_example_gwas_information.txt",header = T)
gwas <- gwas[,c(10,16,19)]
locus <- merge(locus,gwas,by="RSID")
locus$Effect <- abs(locus$Effect)
locus$Correlation <- (locus$Correlation)^2
names(locus)

# Select relevant columns for GRPS calculation
# Columns: gene, Correlation, RSID, peak, Effect, log10P
locus <- locus[,c(8,4,1,18,19,20)]
locus <-locus[!duplicated(locus),]
length(unique(locus$gene))

df <- locus
df <- locus[1,]
df <- df[!duplicated(df),]
set <- df[-1,]
a <- as.character(unique(locus$gene))
info <- locus[1,]
info <- info[!duplicated(info),]
result <- info[-1,]



#Keep the SNP with the max -log10(p-value) in each peak
for (i in 1:length(a)) {
  test <- locus[locus$gene == a[i],]
  head(test)
  peak <- unique(test$peak)
  for (j in 1:length(peak)) {
    ttest <- test[test$peak == peak[j],]
    info <- ttest[which(ttest$log10P == max(ttest$log10P)),]
    info <- info[1,]
    result <- rbind(result,info)
  }
  
}

length(unique(result$gene))
length(unique(locus$gene))

# Calculate pairwise LD using local plink reference panel
library(ieugwasr)
variant <- unique(result$RSID)
test <- ld_matrix_local(variant,
                        bfile = "/Users/zhangjing/Desktop/urate_final/susie/susie/final_merged",
                        plink_bin = "/Users/zhangjing/Desktop/analysis/plink_mac_20250819/plink",
                        with_alleles = F)

# Convert to squared correlation matrix (R²)
ld_mat <- as.matrix(test)
ld_mat <- ld_mat*ld_mat
ld_list <- as.data.frame(as.table(ld_mat))
colnames(ld_list) <- c("SNP_A", "SNP_B", "R2")
ld_list <- ld_list[!duplicated(t(apply(ld_list[, c("SNP_A", "SNP_B")], 1, sort))), ]
ld <- ld_list

# GRPS computation
set <- set[1,]
set <- set[-1,]
locus <- result
colnames(locus)[3] <- "snp"
for (i in 1:length(a)) {
  test <- locus[locus$gene == a[i],]
  head(test)
  
  if (dim(test)[1] == 1) {
    lead <- test[which(test$log10P == max(test$log10P)),]
    lead <-  lead$snp
    ttest <- test
    ttest$module_lead <- lead
    ttest$coef <- 1
    set <- rbind(set,ttest)
    
  } else { 
    lead <- test[which(test$log10P == max(test$log10P)),]
    if (dim(lead)[1] == 1) {
      lead <-  lead$snp
    } else {
      lead <- lead[which(lead$Correlation == max(lead$Correlation)),]
      lead <-  lead$snp
    }
    
    for (j in 1:length(unique(test$snp))) {
      ttest <- test[j,]
      ttest$module_lead <- lead
      snp1 <- ttest$snp
      snp2 <- lead
      coef <- ld[ld$SNP_A == snp1 & ld$SNP_B == snp2 | 
                   ld$SNP_A == snp2 & ld$SNP_B == snp1,]
      if (dim(coef)[1] == 0) {
        ttest$coef <- 1
        set <- rbind(set,ttest)
      } else { 
        ttest$coef <- 1-coef$R2[1]
        set <- rbind(set,ttest)
      }
    }
    
  }
  
}


set$coef[set$snp == set$module_lead] <- 1


final <- as.data.frame(unique(set$gene))
names(final) <- "gene"

head(set)
for (i in 1:length(a)) {
  df <- set[set$gene == a[i],]
  df$GRPS <- df$Correlation*df$Effect*df$coef
  sum <- sum(df$GRPS)
  final[i,2] <- sum
}

names(final) <- c("gene","GRPS")

# Sort genes by GRPS score (descending)
final <- final[order(final$GRPS,decreasing = T),]
head(final)

write.table(final,"./LSS-and-GRPS/output/locus_example_grps_result.txt",quote = F,sep = '\t',row.names = F)



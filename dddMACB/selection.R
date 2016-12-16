
# pass a vcf (ideally from healthy cohort)
# and the set of elements to calculate constraint in
# and variants from unaffected parents

# uses the dddMAPS library

library(stringr)
library(plyr)
library(optparse)
source("../R/annotation_tools.R")
source('../data/mutation_null_model.R')
source("~/software/dddMAPS/dddMAPS/MAPS.R")


### command line options
option_list <- list(
  make_option("--vars", help = "Pass list of variants in unaffected parents"),
  make_option("--elements", help = "Pass elements in which DNMs fall"),
  make_option("--null_model", default = "../data/ddd_synonymous_lm.RData"),
  make_option("--output", help = "location to save data frame with windows stats")
)

args <- parse_args(OptionParser(option_list=option_list))

# load in regions file with required columns: chr, start, stop
elements <- read.table(args$elements, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vars = read.table(args$vars, header = TRUE, sep = "\t")
load(args$null_model)

# add region id
if !("region_id" %in% colnames(elements)){
  elements$region_id = paste0(elements$chr, ":", elements$start, "-", elements$stop)
}


# get sequence context
if !("seq" %in% colnames(elements)) {
  elements$seq = as.character(getSeq(Hsapiens, elements$chr, elements$start, elements$stop))
}


# get probability per element
elements$p_snp_null = 2 * sapply(elements$seq, p_sequence)


# add expected per element
elements$expected = predict(synonymous_lm, elements)


# add observed per element
v = filter_with_bed(vars, elements)
v = subset(v, allele_count < 14)
v$region_id = get_region_id(v, elements)
o = ddply(elements, "region_id", function(df) data.frame(observed = sum(df$allele_count)))
elements$observed = o$observed[match(o$region_id, elements$region_id)]



# add observed/expected
elements$obs_exp_ratio = elements$observed/elements$expected



# calculate Z score from observed and expected
var_Z_score = function(observed, expected){
  Xsq_vals = (observed- expected)^2/expected
  excess = ifelse(observed > expected, -1, 1)
  Z = sqrt(Xsq_vals) * excess
  
  # use trimmed z scores to get standard deviation
  Z_trimmed = Z[ (Z > -5) & (Z < 5)]
  Z_trimmed = Z_trimmed[!is.na(Z_trimmed)]
  Z_sd = sd(Z_trimmed)
  
  # divide ALL Z scores by sd from middle set
  Z_normalized = Z/Z_sd
  
  return(Z_normalized)
}


# add Z score
elements$z_score = var_Z_score(elements$observed, elements$expected)


# add phastcons100 score
library(phastCons100way.UCSC.hg19)
element_intervals = GRanges(seqnames=elements$chr, IRanges(start = elements$start, width = elements$stop - elements$start + 1))
elements$phastcons100 = scores(phastCons100way.UCSC.hg19, element_intervals)


# output annotated elements
write.table(elements, file = output, col.names = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)


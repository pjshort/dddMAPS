# R script to use as a library for doing mutability adjusted proportion of singletons

# important functions for user:
# maps_fit: takes a dataframe of synonymous variants (snps) with columns chr, pos, ref, alt, and allele count, returns maps_lm linear model
# maps_adjust: takes a dataframe of variants with columns chr, pos, ref, alt, allele count, and 'split factor' where split factor is
# a categorization of variants e.g. vep consequence, quantiles of a scoring metric etc.

# used internally: get_context and get_tri

# requires
print("Loading required packages...")
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(plyr)
library(stringr)

print("Loading mutation data and gencode transcripts...")
mu_snp <- read.table("~/phd/code/dddMAPS/data/forSanger_1KG_mutation_rate_table.txt", header=TRUE)

gencode = read.table("~/phd/code/dddMAPS//data/gencode.v19.CDS.probe_overlap.min10_coverage.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
noncoding_intervals = read.table("~/phd/code/de_novo_noncoding/data/de_novo_analysis_regions.noncoding_only.8August2016.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)  # only needed for noncoding analysis
control_introns = read.table("~/phd/code/dddMAPS//data/noncoding_control_elements.10bp_buffer.min10_coverage.30bp_element_minimum.30x_probe_coverage_minimum.no_ddg2p_overlap.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sequences = rbind(gencode[,c("chr", "start", "stop", "seq")], noncoding_intervals[,c("chr", "start", "stop", "seq")], control_introns[,c("chr", "start", "stop", "seq")])

### get the trinucleotide context
 
get_tri = function(interval_idx, pos, intervals) {
  
  start = pos - intervals$start[interval_idx]
  end = pos - intervals$start[interval_idx] + 2

  # explanation:
  # suppose you pick the first base in the sequence, then pos = 0, but R is 1-based
  # so if seq starts at 1000, then 1010 will have pos = 10 (the 10th base) and pos + 2 = 12 (the 12th base)
  # this is the correct tri-nucleotide since 1010 is the 11th base after 1000

  chr = intervals$chr[interval_idx]
  
  seq = as.character(intervals$seq[interval_idx])

  if (start == 0 | end > nchar(seq)) {
    # retrieve using getSeq

    if (!grepl("chr", chr)){
      chr = paste0("chr", chr)
    }
    
    tri = getSeq(Hsapiens, chr, pos - 1, pos + 1)  # uses absolute start/stop


  } else {
    tri = substr(seq, start, end)  # use relative position
  }

}

get_context = function(chromosomes, positions, intervals) {
  # faster version using intervals with $seq column. must be sure that all chr, pos are contained in this set of intervals

  chromosomes = gsub("^chr", "", chromosomes)
  intervals$chr = gsub("^chr", "", intervals$chr)
  
  p = GRanges(seqnames=Rle(chromosomes), ranges = IRanges(start = positions, end = positions))

  if (!is.null(intervals)) {  # the fast way
    i = GRanges(seqnames=Rle(intervals$chr), ranges = IRanges(start = intervals$start, end = intervals$stop))
  } else {
    # the SLOW way using hg19 sequence retrieval for all tri-nucleotides
  }
  
  
  # find overlap between denovos and annotated noncoding
  interval_hits_idx = findOverlaps(p, i, select = "first")
  
  context = mapply(get_tri, interval_hits_idx, positions, MoreArgs = list("intervals" = intervals))
  context = sapply(context, as.character)

  return(context)

}

reverse_complement = function(seq){
  return(as.character(reverseComplement(DNAString(seq))))
}


maps_fit = function(synonymous_vars){
  # synonymous vars should have chr, pos, ref, alt, allele_count
  # mu snps is three columns: from, to, mu_snp
  
  if (!("context" %in% colnames(synonymous_vars))){
    # retrieve context info
    print("Getting tri-nucleotide context for each synonymous variant.")
    
    if (!(any(grepl("chr", synonymous_vars$chr)))) {
      synonymous_vars$chr = paste0("chr", synonymous_vars$chr)
    }
    
    synonymous_vars$context = get_context(synonymous_vars$chr, synonymous_vars$pos, gencode)
  }

  if (!("alt_context" %in% colnames(synonymous_vars))){
    synonymous_vars$alt_context = synonymous_vars$context
    str_sub(synonymous_vars$alt_context, 2, 2) <- synonymous_vars$alt  # replace with alt
  }

  print("Merging the synonymous variants by tri-nucleotide context.")
  synonymous_tri = ddply(synonymous_vars, c('context', 'alt_context', 'ref', 'alt'), function(x) {
    data.frame(n=nrow(x), singletons=sum(x$allele_count == 1), doubletons=sum(x$allele_count == 2),
               tripletons=sum(x$allele_count == 3), quad=sum(x$allele_count == 4), quint=sum(x$allele_count == 5),
               ac_gt_five=sum(x$allele_count > 5), ac_lt_five=sum(x$allele_count < 5), rare_var=sum(x$allele_count < 16))})

  synonymous_tri = merge(synonymous_tri, mu_snp, by.x = c("context", "alt_context"), by.y = c("from", "to"))

  expected_proportion_singleton_lm = lm(singletons/n ~ mu_snp, synonymous_tri, weights = synonymous_tri$n)

  return(expected_proportion_singleton_lm)
}


maps_adjust = function(variants, split_factor, maps_lm, minimum = 100) {
  # take a dataframe of variants that are separated by some identifier called split_factor e.g. "cadd_score"
  # maps_model is a linear model learned on presumed synonymous mutations using lm(singletons/n ~ mu_snp)
  # with mu_snp from the tri-nucleotide mutation model

  reqd_columns = c("chr", "pos", "allele_count")
  if (!all(reqd_columns %in% colnames(variants))){
    warning("One or more required columns missing (chr, pos, allele_count).")
  }
  
  if ("ac" %in% colnames(variants)){
    colnames(variants)[colnames(variants) == "ac"] = "allele_count"
  }
  
  # check if context in colnames. if not, get context for each of the variants and average over them all
  if (!("context" %in% colnames(variants))){
    # retrieve context info
    print("Getting tri-nucleotide context for each input variant.")

    variants$context = get_context(as.character(variants$chr), variants$pos, sequences)
  }

  if (!("alt_context" %in% colnames(variants))){
    variants$alt_context = variants$context
    str_sub(variants$alt_context, 2, 2) <- variants$alt  # replace with alt
  }

  # calculate the singleton_ratio_raw for each split_factor
  variant_split = split(variants, f = split_factor)
  ps_raw = sapply(variant_split, function(d) sum(d$allele_count == 1)/nrow(d))
  counts = sapply(variant_split, function(d) nrow(d))
  standard_error = mapply(function(p,n) sqrt(p*(1-p)/n), ps_raw, counts)

  # calculate average mu_snp for each level in split_factor
  print("Calculating the average mutability across each category.")
  mu_average = sapply(variant_split, function(d) {
    tri = ddply(d, c('context', 'alt_context', 'ref', 'alt'), function(x) {
    data.frame(n=nrow(x), singletons=sum(x$allele_count == 1))})
    tri = merge(tri, mu_snp, by.x = c("context", "alt_context"), by.y = c("from", "to"))
    return(sum(tri$mu_snp*tri$n)/sum(tri$n))})

  print("Predicting the mutability based on sequence context and adjusting...")

  ps_predicted = predict(maps_lm, data.frame(mu_snp = unlist(mu_average)))
  ps_adjusted = ps_raw - ps_predicted
  
  print("Removing any data points below minimum... Default is 100.")
  ps_adjusted[counts < minimum] = NA
  standard_error[counts < minimum] = NA
  
  return(list("ps_adjusted" = ps_adjusted, "standard_error" = standard_error))
}

maps_grid <- function(variants, split_factor_x, split_factor_y, maps_lm){
  # loop over the split_factor_y - this gives us a maps_adjust output for each level of y
  variants_by_y = split(variants, variants[,split_factor_y])
  
  ps_mat = matrix(ncol = length(levels(variants[,split_factor_x])), nrow = length(levels(variants[,split_factor_y])))
  colnames(ps_mat) = levels(variants[,split_factor_x])
  rownames(ps_mat) = levels(variants[,split_factor_y])
  
  for (i in seq_along(variants_by_y)) {
    print(i)
    v = variants_by_y[[i]]
    ps_mat[i,] = maps_adjust(v, v[,split_factor_x], maps_lm, minimum = 100)$ps_adjusted
  }
}

ps_raw = function(variants, split_factor){
  
  # calculate the singleton_ratio_raw for each split_factor
  variant_split = split(variants, f = split_factor)
  ps_raw = sapply(variant_split, function(d) sum(d$allele_count == 1)/nrow(d))
  counts = sapply(variant_split, function(d) nrow(d))
  standard_error = mapply(function(p,n) sqrt(p*(1-p)/n), ps_raw, counts)
  
  return(list("ps_raw" = ps_raw, "standard_error" = standard_error))
}

maps_ggplot = function(split_levels, ps_adjusted, standard_error, already_ordered = FALSE, add_coding_fixed = c(0, 0.0548, 0.141), add_synonymous_fixed = FALSE, colors = NULL, score_name = "Scoring Metric"){
  # makes a simple ggplot of the mutability adjusted prop of singletons with error bars

  
  df = data.frame(split_level = split_levels, ratio = ps_adjusted, se = standard_error)
  
  if (!is.null(colors)){
    df$colors = colors
  } else {
    df$colors = NULL
  }
  
  if (!already_ordered){
    df = df[order(df$ratio),]
  }
  
  if (any(add_coding_fixed != FALSE)){  # these are from VEP consequences on DDD unaffected parents
    coding_fixed = data.frame(split_level = c("Synonymous", "Missense", "Stop Gained"), ratio = add_coding_fixed, se = c(0.0007, 0.0006, 0.0035))
    coding_fixed$colors = "Coding VEP"
    df$colors = score_name
    df = rbind(df, coding_fixed)
    df$colors = factor(df$colors, levels = c(score_name, "Coding VEP"), ordered = TRUE)
    colors = TRUE
  } else if (add_synonymous_fixed == TRUE) {
    coding_fixed = data.frame(split_level = "Synonymous", ratio = 0, se = 0.0007)
    coding_fixed$colors = "Synonymous"
    df$colors = score_name
    df = rbind(df, coding_fixed)
    colors = TRUE
  }

  print("Removing NaNs (insufficient counts).")
  df = df[!is.nan(df$ratio),]

  df$split_level = factor(df$split_level, levels = df$split_level)
  
  if (is.null(colors)){  # plot all in black
    limits = aes(ymin = df$ratio - 1.96*df$se, ymax = df$ratio + 1.96*df$se)
    ggplot(df, aes(split_level, ratio)) +
      geom_pointrange(limits, size = 1.25) + coord_flip() +
      xlab(score_name) + ylab("Mutability Adjusted Proportion of Singletons") +
      theme_bw(base_size = 18) + 
      theme(strip.text = element_text(color="black"),strip.background = element_rect(fill="white", size=0),panel.border = element_blank()) + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
      theme(legend.title = element_blank())
  } else {  # plot with colored entries
    limits = aes(ymin = df$ratio - 1.96*df$se, ymax = df$ratio + 1.96*df$se)
    ggplot(df, aes(split_level, ratio, color = colors)) +
      geom_pointrange(limits, size = 1.25) + coord_flip() +
      xlab(score_name) + ylab("Mutability Adjusted Proportion of Singletons") +
      theme_bw(base_size = 18) + 
      theme(strip.text = element_text(color="black"),strip.background = element_rect(fill="white", size=0),panel.border = element_blank()) + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
      theme(legend.title = element_blank())
  }
}


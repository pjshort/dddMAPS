library(plyr)
library(stringr)

# R script to use as a library for doing mutability adjusted proportion of singletons

# important functions for user:
# maps_fit: takes a dataframe of synonymous variants (snps) with columns chr, pos, ref, alt, and allele count, returns maps_lm linear model
# maps_adjust: takes a dataframe of variants with columns chr, pos, ref, alt, allele count, and 'split factor' where split factor is
# a categorization of variants e.g. vep consequence, quantiles of a scoring metric etc.

# used internally: get_context and get_tri

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

mu_snp <- read.table("../data/forSanger_1KG_mutation_rate_table.txt", header=TRUE)
gencode = read.table("../data/gencode_protein_coding_genes_v19_+strand.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
CNEs = read.table("../data/noncoding_regions.txt", header = TRUE, sep = "\t")


### get the trinucleotide context

get_tri = function(interval_idx, pos, intervals) {
  start = pos - intervals$start[interval_idx]
  end = pos - intervals$start[interval_idx] + 2

  # explanation:
  # suppose you pick the first base in the sequence, then pos = 0, but R is 1-based
  # so if seq starts at 1000, then 1010 will have pos = 10 (the 10th base) and pos + 2 = 12 (the 12th base)
  # this is the correct tri-nucleotide since 1010 is the 11th base after 1000

  chr = intervals$chr[interval_idx]
  if (!grepl("chr", chr)){
    chr = paste0("chr", chr)
  }
  seq = as.character(intervals$seq[interval_idx])

  if (start == 0 | end > nchar(seq)) {
    # retrieve using getSeq

    tri = getSeq(Hsapiens, chr, pos - 1, pos + 1)  # uses absolute start/stop


  } else {
    tri = substr(seq, start, end)  # use relative position
  }

}

get_context = function(chromosomes, positions, intervals) {
  # faster version using intervals with $seq column. must be sure that all chr, pos are contained in this set of intervals

  p = GRanges(seqnames=Rle(chromosomes), ranges = IRanges(start = positions, end = positions))

  if (!is.null(intervals)) {  # the fast way
    i = GRanges(seqnames=Rle(intervals$chr), ranges = IRanges(start = intervals$start, end = intervals$stop))
  } else {
    # the SLOW way using hg19 sequence retrieval for all tri-nucleotides
  }

  # find overlap between denovos and annotated CNEs
  interval_hits_idx = findOverlaps(p, i, select = "first")

  context = mapply(get_tri, interval_hits_idx, positions, MoreArgs = list("intervals" = intervals))
  context = sapply(context, as.character)

  return(context)

}


maps_fit = function(synonymous_vars){
  # synonymous vars should have chr, pos, ref, alt, allele_count
  # mu snps is three columns: from, to, mu_snp

  if (!("context" %in% colnames(synonymous_vars))){
    # retrieve context info
    print("Getting tri-nucleotide context for each synonymous variant.")
    synonymous_vars$context = get_context(paste0("chr", synonymous_vars$chr), synonymous_vars$pos, gencode)
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


maps_adjust = function(variants, split_factor, maps_lm, noncoding = TRUE) {
  # take a dataframe of variants that are separated by some identifier called split_factor e.g. "cadd_score"
  # maps_model is a linear model learned on presumed synonymous mutations using lm(singletons/n ~ mu_snp)
  # with mu_snp from the tri-nucleotide mutation model


  # check if context in colnames. if not, get context for each of the variants and average over them all
  if (!("context" %in% colnames(variants))){
    # retrieve context info
    print("Getting tri-nucleotide context for each input variant.")
    if (noncoding == TRUE){
      variants$context = get_context(as.character(variants$chr), variants$pos, CNEs)
    } else {
      variants$context = get_context(paste0("chr", variants$chr), variants$pos, gencode)
    }
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

  return(list("ps_adjusted" = ps_adjusted, "standard_error" = standard_error))
}



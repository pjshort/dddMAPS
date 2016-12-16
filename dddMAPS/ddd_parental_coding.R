# example of MAPS.R usage on DDD coding variants
source("./MAPS.R")
source("./MAPS_plotting_extras.R")
library(Biostrings)
library(stringr)

# run with pre-retrieved context
parental_gencode_snps = read.table("../data/unaffected_parent_alleles_all_chromosomes.CADD.gencode_v19_protein_coding.txt", header = TRUE, sep = "\t")

parental_gencode_snps$alt_context = parental_gencode_snps$context
str_sub(parental_gencode_snps$alt_context, 2, 2) <- parental_gencode_snps$alt  # replace with alt

trinucleotide_rev_complements = read.table("../data/trinucleotide_mapping.txt", header = FALSE, sep = "\t")

mutation = paste0(parental_gencode_snps$context, ">", parental_gencode_snps$alt_context)
mutation[mutation %in% trinucleotide_rev_complements$V2] = as.character(trinucleotide_rev_complements$V1)[!is.na(match(mutation,trinucleotide_rev_complements$V2))]

synonymous_parental_vars = subset(parental_gencode_snps, vep_consequence == "synonymous_variant")

# fit linear model to synonymous variants
maps_lm = maps_fit(synonymous_parental_vars)

# split cadd score into bins
parental_gencode_snps$cadd_bin = cut(parental_gencode_snps$scaled_CADD, c(seq(0,40,2), 50))

# get mutability adjusted proportion of singletons 
out = maps_adjust(variants = parental_gencode_snps, split_factor = parental_gencode_snps$cadd_bin, maps_lm = maps_lm)

# plot MAPS for each cadd bin
maps_ggplot(names(out$ps_adjusted), out$ps_adjusted, out$standard_error, already_ordered = FALSE, score_name = "CADD")

# add fixed coding:
maps_ggplot(names(out$ps_adjusted), out$ps_adjusted, out$standard_error, already_ordered = FALSE, add_coding_fixed = TRUE)

# add bar plot of counts
maps_plus_bar(names(out$ps_adjusted), out$ps_adjusted, out$standard_error, parental_gencode_snps, parental_gencode_snps$cadd_bin, add_coding_fixed = TRUE)


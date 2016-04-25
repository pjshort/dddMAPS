# example of MAPS.R usage on DDD coding variants
source("./MAPS.R")
library(Biostrings)
library(stringr)

# run with pre-retrieved context
parental_gencode_snps = read.table("../data/gencode_parental_snps_+context.txt", header = TRUE, sep = "\t")
parental_gencode_snps$alt_context = parental_gencode_snps$context
str_sub(parental_gencode_snps$alt_context, 2, 2) <- parental_gencode_snps$alt  # replace with alt

trinucleotide_rev_complements = read.table("../data/trinucleotide_mapping.txt", header = FALSE, sep = "\t")

mutation = paste0(parental_gencode_snps$context, ">", parental_gencode_snps$alt_context)
mutation[mutation %in% trinucleotide_rev_complements$V2] = as.character(trinucleotide_rev_complements$V1)[!is.na(match(mutation,trinucleotide_rev_complements$V2))]

#run these lines instead of above to generate the context (takes longer)
#parental_gencode = read.table("../data/gencode_parental_alleles_FULL.txt", header = TRUE, sep = "\t")
#parental_gencode = subset(parental_gencode, nchar(as.character(ref)) == 1 & nchar(as.character(alt)) == 1)

synonymous_parental_vars = subset(parental_gencode_snps, vep_consequence == "synonymous_variant")


# fit linear model to synonymous variants
maps_lm = maps_fit(synonymous_parental_vars)

# get mutability adjusted proportion of singletons for all coding variants
out = maps_adjust(variants = parental_gencode_snps, split_factor = parental_gencode_snps$vep_consequence, maps_lm = maps_lm, noncoding = FALSE)

# plot MAPS for each vep consequence
maps_ggplot(names(out$ps_adjusted), out$ps_adjusted, out$standard_error, already_ordered = FALSE)

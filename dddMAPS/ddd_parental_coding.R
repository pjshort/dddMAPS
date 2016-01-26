# example of MAPS.R usage on DDD coding variants

source("./MAPS.R")

parental_gencode = read.table("../data/gencode_parental_alleles_FULL.txt", header = TRUE, sep = "\t")
parental_gencode = subset(parental_gencode, nchar(as.character(ref)) == 1 & nchar(as.character(alt)) == 1)
synonymous_parental_vars = subset(parental_gencode, vep_consequence == "synonymous_variant")

# fit linear model to synonymous variants
maps_lm = maps_fit(synonymous_parental_vars)

# get mutability adjusted proportion of singletons for all coding variants
out = maps_adjust(variants = parental_gencode, split_factor = parental_gencode$vep_consequence, maps_lm = maps_lm, noncoding = FALSE)

# plot MAPS for each vep consequence
ggplot_MAPS(names(out$ps_adjusted), out$ps_adjusted, out$standard_error, already_ordered = FALSE)
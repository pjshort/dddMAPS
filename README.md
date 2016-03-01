# dddMAPS
Mutability adjusted proportion of singletons (MAPS) calculation using synonymous SNPs from the DDD unaffected parents.

This method was used in the ExAC 60k genomes paper and has been adopted here for use in analyzing selection in coding/non-coding regions for the DDD cohort.

# Requirements

GenomicRanges
BSgenome.Hsapiens.UCSC.hg19 (not strictly essential - download will take a few minutes.)

```R
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
```

The code will run much faster with the Ensembl transcripts file that has pre-loaded sequence context for each transcript.

# Fit maps_fit and run maps_adjust

There are two main relevant functions, maps_fit and maps_adjust. A maps_lm (MAPS linear model) has already been fit to the DDD unaffected parents and is included in the data folder as ddd_parents_maps_lm.RData. Otherwise, a maps_lm can be fit with:

```R
synonymous_vars = # dataframe with chr, pos, ref, alt, allele_count in a presumed healthy population
maps_lm = maps_fit(synonymous_vars)
```

The maps_fit function will look for a 'context' column in synonymous vars. If it is not found, it will look for the gencode_protein_coding_genes_v19_+strand.txt file in the data directory to quickly retrieve the context. If this is not found, it will take the long way, using BSgenome.Hsapiens.UCSC.hg19 to get the context.

Mutability adjusted proportion of singletons can be calculated for a set of variants in a dataframe organized/separated into groups with the same required columns as synonymous variants (chr, pos, ref, alt, allele_count) and an additional column specifying a categorization or 'split_factor' using maps_adjust.

For instance, if a dataframe of protein-coding variants are annotated with a vep_consequence:

```R
split_factor = variants$vep_consequence
out = maps_adjust(variants, split_factor, maps_lm, noncoding = FALSE)
ps_adjusted = out$ps_adjusted
standard_error = out$standard_error
```

ps_adjusted will be a named vector with the mutability adjusted proportion of singletons for each vep_consequence. If you include synonymous variants, you will see the maps value will be very close to zero which is expected (nearly all singleton variants in these regions are fit by underlying variation in mutability).

# Plotting Results

Finally, a simple maps_ggplot function is included to sort and plot the MAPS values by split_factor:
```R
maps_ggplot = function(names(ps_adjusted), ps_adjusted, standard_error, already_ordered = FALSE)
```

# Coding Region Example
See dddMAPS/ddd_parental_coding.R for top-to-bottom example for the DDD parental coding regions.

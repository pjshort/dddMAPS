The dddMAPS project contains two libraries that can be used to calculate constraint in coding/noncoding regions. The first, dddMAPS, is focused on 'mutability adjusted proportion of singletons', a measure of purifying selection. The second, dddMACB, is focused on measuring constraint by comparing the observed amount of damaging variation in an element (as scored by CADD) compared to the amount expected under a selection-neutral model.

# dddMAPS
Mutability adjusted proportion of singletons (MAPS) calculation using synonymous SNPs from the DDD unaffected parents.

This method was used in the ExAC 60k genomes paper and has been adopted here for use in analyzing selection in coding/non-coding regions for the DDD cohort.

## Requirements

GenomicRanges
BSgenome.Hsapiens.UCSC.hg19 (not strictly essential - download will take a few minutes.)

```R
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
```

The code will run much faster with the Ensembl transcripts file that has pre-loaded sequence context for each transcript. This is on the sanger farm - please get in touch for the path to these files!

## Fit maps_fit and run maps_adjust

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

## Plotting Results

Finally, a simple maps_ggplot function is included to sort and plot the MAPS values by split_factor:
```R
maps_ggplot = function(names(ps_adjusted), ps_adjusted, standard_error, already_ordered = FALSE)
```

## Coding Region Example
See dddMAPS/ddd_parental_coding.R for top-to-bottom example for the DDD parental coding regions.


#dddMACB
This stands for 'mutability adjusted cadd burden' and is being used as a simple and flexible way to calculate constraint scores that are comparable between coding and non-coding regions.

The workflow for calculating MACB is simply:
1. Generate an 'exhaustive allele file' for all regions of interest. This includes every ref/alt combo at each position across the region.
2. Calculate CADD scores for each of these positions.
3. Use the mutability and the CADD score to calculate the weighted average of CADD scores across the region/element provided. This can be thought of as the expected value of the CADD score for a de novo mutation in the region.
4. 

## Generate Exhaustive SNP file

Regions file requires chr, start, stop, and 'region_id'. This could be the name of the gene, or in the case of noncoding elements it is simply 'chr:start-stop'. This is used to calculate the expected mutability for a given region. If a region id is not provided, the script will look for 'gene' column and if not found, a region_id column will be created.

```bash
bsub -J exhaustive_noncoding[1-7540:10] -q normal -R'select[mem>300] rusage[mem=300]' -M300 \
-o $pjs/MACB/logs/non_coding_exhaustive.%I.out /software/R-3.2.2/bin/Rscript create_exhaustive_allele_files.Rscript \ --index=\$LSB_JOBINDEX --index_step=10 --regions=~/reference_data/CNEs_subtract_CDS.txt \ --out_base=$pjs/MACB/alleles/noncoding_exhaustive_allele
```

## Calculate CADD scores
This requires the CADD SNP score tabix file.

```bash
bsub -q normal -J "noncoding_exhaustive_cadd[1-7540:10]" -R'select[mem>100] rusage[mem=100]' -M100 -o \ /lustre/scratch113/projects/ddd/users/ps14/CADD/logs/noncoding_exhaustive.%I.out python -u \ ~/software/SingletonMetric/python/TabixScores.py --tabix /lustre/scratch113/projects/ddd/users/ps14/CADD/whole_genome_SNVs.tsv.gz \ --variants $pjs/MACB/alleles/noncoding_exhaustive_allele.\$LSB_JOBINDEX.txt \
--variants_out $pjs/MACB/alleles/noncoding_exhaustive_allele.$LSB_JOBINDEX.CADD.txt \
--score CADD
```

## Calculate MACB for each of the regions.

```bash
bsub -J MACB_noncoding -q normal -R'select[mem>10000] rusage[mem=10000]' -M10000 -o \
$pjs/MACB/logs/noncoding_MACB.%I.out /software/R-3.2.2/bin/Rscript \
~/software/dddMAPS/dddMACB/calculate_MACB_null.Rscript \
--input_base=$pjs/MACB/alleles/noncoding_alleles_exhaustive \
--index_start=1 --index_stop=7540 --index_step=10 --out=$pjs/MACB/non_coding_elements_MACB.txt \
--score=MACB
```

To run a MACB that is more like loss-of-function intolerance, MACB25 calculates the expected proportion of sites with CADD >25 weighted by mutability. By multiplying this by the estimate of rare variants per element, this provides the expected number of LoF-like variants in an element.


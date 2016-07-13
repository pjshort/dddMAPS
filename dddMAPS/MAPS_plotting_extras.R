# plotting supplements for MAPS.R

library(gridExtra)
library(ggplot2)
library(plyr)
library(stringr)

maps_plus_bar = function(split_levels, ps_adjusted, standard_error, variants, split_factor, add_coding_fixed = FALSE, add_synonymous_fixed = FALSE, colors = NULL, score_name = "Scoring Metric"){
  # makes mutability adjusted proportion of singletons plot by split_levels plus bar plot with number of variants per bin
  
  # get the main singleton plot
  main = maps_ggplot(split_levels, ps_adjusted, standard_error, add_coding_fixed = add_coding_fixed, add_synonymous_fixed = add_synonymous_fixed, already_ordered = TRUE, colors = colors, score_name = score_name)
  
  # get the number of variants per split level
  counts = sapply(split(variants, split_factor), nrow)
  counts = data.frame(counts = counts, split_level = split_levels)
  
  if (add_coding_fixed == TRUE){  # these are from VEP consequences on DDD unaffected parents
    # coding_counts = data.frame(split_level = c("Synonymous", "Missense", "Stop Gained"), counts = c(486901,870066,22830))
    coding_counts = data.frame(split_level = c("Synonymous", "Missense", "Stop Gained"), counts = c(0,0,0))
    counts = rbind(counts, coding_counts)
    counts$split_level = factor(counts$split_level, levels = c(as.character(counts$split_level), "Synonymous", "Missense", "Stop Gained"))
  } else if (add_synonymous_fixed == TRUE) {
    coding_counts = data.frame(split_level = "Synonymous", counts = 0)
    counts = rbind(counts, coding_counts)
    counts$split_level = factor(counts$split_level, levels = c(as.character(counts$split_level), "Synonymous"))
  } else {
    counts$split_level = factor(counts$split_level, levels = c(split_levels))
  }
  
  counts = ggplot(counts) + geom_bar(aes(split_level, counts/1000), stat = "identity") + coord_flip() + xlab("Number of Variants (Thousands)") + ylab("Counts") +
    theme_bw(base_size = 18) + 
    theme(strip.text = element_text(color="black"),strip.background = element_rect(fill="white", size=0),panel.border = element_blank()) + 
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    theme(legend.title = element_blank())
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  grid.arrange(counts, main, ncol = 2, widths=c(2, 3.5))
  
}
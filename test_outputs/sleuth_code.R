library(sleuth)
library(dplyr)
#set working directory
setwd('~/miniProject/test_outputs')
#read csv file into stab
stab <- read.csv('input_table.csv', header = TRUE, stringsAsFactors = FALSE)
#initialize sleuth object
so <- sleuth_prep(stab)

#===== Differential expression analysis (2dpi vs 6dpi) =====

#fit model to compare conditions
so <- sleuth_fit(so, ~condition, 'full')
#fit reduced model to compare in likelihood ratio test
so <- sleuth_fit(so, ~1, 'reduced') 
#perform likelihood ratio test for diffy exp. b/w conditions
so <- sleuth_lrt(so, 'reduced', 'full') 


#===== Extract & find sig. results =====
sleuth_tab <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
#filter FDR<0.05 transcripts & sort by pval
sleuth_sig <- dplyr::filter(sleuth_tab, qval <= 0.05) %>% dplyr::arrange(pval) 
#select columns want to output
sleuth_sig <- dplyr::select(sleuth_sig, target_id, test_stat, pval, qval)
#write tab-delimited table to log file
write.table(sleuth_sig, file = 'miniProject.log', quote = FALSE, row.names = FALSE, sep = '\t',append = TRUE)


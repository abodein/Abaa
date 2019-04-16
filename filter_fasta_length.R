args <- commandArgs(TRUE)

library(tidyverse)
library(ShortRead)

fasta <- readFasta(args[1])

length_freq <- width(fasta) %>% table %>% as.data.frame() %>% set_names("length", "freq") %>% mutate(length =length %>% as.character() %>% as.numeric())
good_seq <- fasta[width(fasta) %in% (length_freq %>% filter(freq >= sd(freq)) %>% pull(length))]

writeFasta(object = good_seq, file = args[2], mode = "w")

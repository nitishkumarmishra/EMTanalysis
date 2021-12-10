library(GenomicFeatures)
library(dplyr)


ensemble_gtf <- makeTxDbFromGFF('../../RSEM_COUNTS_OLD/Counts/Mus_musculus.GRCm39.104.gtf.gz', organism = "Homo sapiens")
fiveUTRs          <- fiveUTRsByTranscript(ensemble_gtf,use.names=TRUE)
length_fiveUTRs   <- width(ranges(fiveUTRs))
the_lengths        <- as.data.frame(length_fiveUTRs)
the_lengths        <- the_lengths %>% group_by(group, group_name) %>% summarise(sum(value))
the_lengths        <- unique(the_lengths[,c("group_name", "sum(value)")])
colnames(the_lengths) <- c("ENST", "5' UTR Length")
the_lengths <- as.data.frame(the_lengths)



library('rtracklayer')
GTF="C:/Users/nmishra/Dropbox/PC/Desktop/Total RNAseq EMT/RSEM_COUNTS_OLD/Counts/Mus_musculus.GRCm39.104.gtf.gz"
GTF <- import(GTF)
GRCm39.104.gtf <- as.data.frame(GTF)

GRCm39.104.selected <- GRCm39.104.gtf %>%
  filter(type=="gene") %>%
  rename(Chr=seqnames, Type=type,ENSG=gene_id, Symbol=gene_name, GeneType=gene_biotype, GeneLength=width) %>%
  select(Chr, Type, GeneType, ENSG, Symbol, GeneLength) 

GRCm39.104.ENST <- GRCm39.104.gtf %>%
  rename(ENSG=gene_id, ENST=transcript_id) %>%
  select(ENSG, ENST) %>%
  tidyr::drop_na()%>%  #filter(!is.na(ENST)) ##Both filter/is.na will work
  distinct()


UTR_Length <- merge(GRCm39.104.ENST, the_lengths, by ="ENST")

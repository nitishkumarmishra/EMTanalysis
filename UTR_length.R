library(GenomicFeatures)
library(dplyr)

refSeq             <- makeTxDbFromUCSC(genom="mm39",tablename="refGene")                     
fiveUTRs          <- fiveUTRsByTranscript(refSeq, use.names=TRUE)
length_fiveUTRs   <- width(ranges(fiveUTRs))
the_lengths        <- as.data.frame(length_fiveUTRs)
the_lengths        <- the_lengths %>% group_by(group, group_name) %>% summarise(sum(value))
the_lengths        <- unique(the_lengths[,c("group_name", "sum(value)")])
colnames(the_lengths) <- c("RefSeq Transcript", "5' UTR Length")


ensemble <- makeTxDbFromEnsembl(organism="Mus musculus",
                                release=NA,
                                circ_seqs=NULL,
                                server="ensembldb.ensembl.org",
                                username="anonymous", password=NULL, port=0L,
                                tx_attrib=NULL)
fiveUTRs          <- fiveUTRsByTranscript(ensemble,use.names=TRUE)
length_fiveUTRs   <- width(ranges(fiveUTRs))
the_lengths        <- as.data.frame(length_fiveUTRs)
the_lengths        <- the_lengths %>% group_by(group, group_name) %>% summarise(sum(value))
the_lengths        <- unique(the_lengths[,c("group_name", "sum(value)")])
colnames(the_lengths) <- c("ensemble Transcript", "5' UTR Length")




GENCODE_gtf <- makeTxDbFromGFF('../../gencode.vM27.annotation.gtf', organism = "Homo sapiens")
fiveUTRs          <- fiveUTRsByTranscript(GENCODE_gtf,use.names=TRUE)
length_fiveUTRs   <- width(ranges(fiveUTRs))
the_lengths        <- as.data.frame(length_fiveUTRs)
the_lengths        <- the_lengths %>% group_by(group, group_name) %>% summarise(sum(value))
the_lengths        <- unique(the_lengths[,c("group_name", "sum(value)")])
colnames(the_lengths) <- c("ensemble Transcript", "5' UTR Length")


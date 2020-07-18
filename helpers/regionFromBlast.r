suppressMessages(library("GenomicRanges"))

args = commandArgs(trailingOnly=TRUE)

blast = args[1]
out = args[2]

min_blast_length = 100

# load blast results
blast_result <- read.delim(blast, head = F)
colnames(blast_result) <- c("qseqid", 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen','qcovs')
blast_result <- subset(blast_result, blast_result$length > min_blast_length)
region_length <- blast_result$qlen[1]

# construct granges
subject_df <- data.frame(seqnames=blast_result$sseqid, start=blast_result$sstart, end = blast_result$send)
subject_mat <- apply(subject_df, 1, function(x) {start=min(x[2:3]); end=max(x[2:3]); x[2] = start; x[3] = end; return(x)})
subject_df <- as.data.frame(t(subject_mat))
subject_gr <- GRanges(subject_df)
subject_gr <- reduce(subject_gr)

# find contig with most massive region
subject_contigs <- levels(seqnames(subject_gr))
mass <- list()
for(contig in subject_contigs){
  subject_contig_gr <- subset(subject_gr, seqnames(subject_gr)==contig)
  mass[[contig]] <- sum(width(subject_contig_gr))  
}
the_contig <- names(mass)[which.max(unlist(mass))]

# find region centr of mass
subject_contig_gr <- subset(subject_gr, seqnames(subject_gr)==the_contig)
subject_contig_df <- data.frame(subject_contig_gr)
subject_contig_df$mid <- (subject_contig_df$start+subject_contig_df$end)/2
M <- sum(subject_contig_df$width)
centr_mass <- (sum(subject_contig_df$mid*subject_contig_df$width))/M
the_region <- GRanges(the_contig, IRanges(start=(centr_mass-region_length), end=(start=centr_mass+region_length))) # mid +- region_length
write.table(data.frame(the_region), quote = F, row.names = F, col.names = F, file = out, sep="\t")

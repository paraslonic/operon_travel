library("data.table")
library("IRanges")

args = commandArgs(trailingOnly=TRUE)
operon = args[1]

## load data
faa_path <- paste0("orthosnake/",operon,"/faa/")
ref_genome <- list.files(faa_path,"*.fasta",full.names = F)[1]
genes_table <- fread(cmd=paste0("grep '>' ",paste0(faa_path,ref_genome)))
colnames(genes_table) <- c("genome","id","product","contig","start","end")

ref_genome_name <- gsub(".fasta", "", ref_genome)
regions <- read.delim(paste0("blast/",operon,"/", ref_genome_name,".blast"), head = F)
colnames(regions) <- c("qseqid", 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

## Find genes
overlapped_genes_index_total <- c()

for(i in 1:nrow(regions)){
  region <- regions[i,]
  genes_table_subset <- subset(genes_table, genes_table$contig == region$sseqid)  
  region_ir <- IRanges(start = region$sstart, end = region$send)
  genes_ir <- IRanges(start = genes_table_subset$start, end = genes_table_subset$end)
  overlaps <- findOverlaps(genes_ir, region_ir, type="within")
  overlapped_genes_index <- as.data.frame(overlaps)[,1]
  overlapped_genes_index_total <- c(overlapped_genes_index_total, overlapped_genes_index)
}

genes_in_region <- apply(genes_table[overlapped_genes_index_total,], 1,  paste,collapse="|")
genes_in_region <- gsub("^>","",genes_in_region)

## Find OGs
og_table <- fread(paste0("orthosnake/",operon,"/Results/Orthogroups.tsv"))

ref_col <- which(colnames(og_table) == ref_genome_name)

og_in_region_index <- which(unlist(og_table[, ref_col, with=F]) %in% genes_in_region)

og_in_region <- og_table[og_in_region_index, 1, with=F]
write.table(og_in_region, paste0("tmp/og_in_region/",operon,"/og_list"), row.names = F, col.names = F, quote = F)


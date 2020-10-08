set -eo pipefail
folder=$1
mkdir -p $folder/tmp 
cat $folder/ffn/*.fasta > $folder/tmp/all_genes_nuc.fasta
cat $folder/faa/*.fasta > $folder/tmp/all_genes_prot.fasta
mkdir -p $folder/Results/ortho/
mkdir -p $folder/Results/ortho/coreogs_nuc
mkdir -p $folder/Results/ortho/coreogs_prot

perl scripts/splitToOg.pl $folder/Results/Orthogroups.txt $folder/tmp/all_genes_nuc.fasta $folder/Results/ortho/coreogs_nuc $folder/Results/Orthogroups_SingleCopyOrthologues.txt >> $folder/plog

perl scripts/splitToOg.pl $folder/Results/Orthogroups.txt $folder/tmp/all_genes_prot.fasta $folder/Results/ortho/coreogs_prot $folder/Results/Orthogroups_SingleCopyOrthologues.txt >> $folder/plog
     
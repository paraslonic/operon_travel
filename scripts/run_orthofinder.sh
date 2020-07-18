folder=$1
threads=$2
orthofinder -t $threads -a $threads -og -f $folder/faa
mkdir -p tmp
find $folder/faa -name 'Orthogroups.txt' -exec cp {} $folder/Results \; 
find $folder/faa -name 'Orthogroups_SingleCopyOrthologues.txt' -exec cp {} $folder/Results \;

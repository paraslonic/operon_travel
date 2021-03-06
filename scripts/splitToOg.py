import sys
import os
from Bio import SeqIO

orthologousGroupsFile = sys.argv[1] 
sequencesFile= sys.argv[2];
ogListFile=sys.argv[3]
outputFolder=sys.argv[4]

try:
    os.makedirs(outputFolder)
except OSError:
    print ("Creation of the directory %s failed" % outputFolder)
else:
    print ("Successfully created the directory %s" % outputFolder)

ogList = []

with open(ogListFile, "r") as oglistfile:
    for og in oglistfile: 
        ogList.append(og.strip())

gene_og = {}

with open(orthologousGroupsFile, "r") as ogfile:
    for ogline in ogfile:
        ogline = ogline.strip()
        [og, genesStr] = ogline.split(": ")
        genes = genesStr.split(" ")
        if not og in ogList: continue
        for gene in genes:
            genome_geneid = tuple(gene.split("|")[0:2])
            gene_og[genome_geneid] = og        

og_records = {}

for seq_record in SeqIO.parse(sequencesFile, "fasta"):
    seq_id = seq_record.id
    seq_genome_geneid = tuple(seq_id.split("|")[0:2])
    if seq_genome_geneid in gene_og:
        og = gene_og[seq_genome_geneid]
        if(og in og_records):
            og_records[og].append(seq_record)
        else:
            og_records[og] = []

for og in og_records:
    fname=outputFolder+"/"+og+".fasta"
    with open(fname,"w") as outfile:
        SeqIO.write(og_records[og],outfile,"fasta")
import os
import re
import glob

results = []

genomes,=glob_wildcards("fna/{genome}.fna")
#genomes = genomes[0:60]

#results.append(expand("coverage/pdu/{genome}.cov", genome=genomes))
#results.append(expand("blast/{genome}.blast", genome=genomes))

#results.append("genomes_with_regions/pdu")
results.append("tmp/result/pdu/out")
results.append("orthosnake/pdu/Results/coreogs_nucleotide.treefile")
results.append("tmp/og_in_region/pdu/og_list")

print("Genomes: ", len(genomes))

rule all:
    input: results

rule blast:
    input: region="regions/{region}.fasta",genome="fna/{genome}.fna"
    output: "blast/{region}/{genome}.blast"
    conda: "envs/env.yaml"
    shell: 
        """
        blastn -query {input.region} -subject {input.genome} -outfmt '6 std qlen qcovs' -out {output}
        """

rule calculate_coverage:
    input: "blast/{region}/{genome}.blast"
    output: "coverage/{region}/{genome}.cov"
    shell: 
        """
        printf "{wildcards.genome}\t" > {output}
        if [ -s "{input}" ] 
        then
            head -1  {input} | awk '{{print $NF}}'  >> {output} # assume qcovs is in last column
        else
            printf "0" >> {output}
        fi
        """

checkpoint select_genomes_with_regions:
    input:  lambda wildcards: 
            ["coverage/{0}/{1}.cov".format(wildcards.region, genome) for genome in genomes]
    output: outdir=directory("genomes_with_regions/{region}")
    params: min_coverage=90
    run:
        region = wildcards.region
        os.mkdir(output.outdir)
        for cov_file in input:
            if(os.path.getsize(cov_file) < 1): continue
            with open(cov_file, "r") as infile:
                line = infile.readline()
                print(line)
                [genome, coverage] = re.split(r'\t+', line)
                if(int(coverage) > params.min_coverage): 
                    os.symlink("../../fna/"+genome+".fna", "genomes_with_regions/"+region+"/"+genome+".fna")


#============================================ ORTHO SNAKE --

def get_selected_genomes(wildcards):
    selected_genomes_folder = checkpoints.select_genomes_with_regions.get(**wildcards).output[0]
    genomes = glob_wildcards(os.path.join(selected_genomes_folder, "{genome}.fna")).genome
    return genomes


rule check_fna:
    input: "genomes_with_regions/{region}/{genome}.fna"
    output: "orthosnake/{region}/fna_final/{genome}.fna"
    conda: "envs/scripts.yaml"
    shell: "python scripts/cropHeader.py -input  {input} -out {output} -n 20"
    
rule prokka:
    input: "orthosnake/{region}/fna_final/{genome}.fna"
    output: dir=directory("orthosnake/{region}/prokka/{genome}"),
            gbk="orthosnake/{region}/prokka/{genome}/{genome}.gbk"
    threads: 4
    conda: "envs/prokka.yaml"
    shell:
        """
        prokka --cpus {threads} --outdir {output.dir} --force --prefix {wildcards.genome} --locustag {wildcards.genome} {input} 2>/dev/null
        """

rule make_faa:
    input: "orthosnake/{region}/prokka/{genome}/{genome}.gbk"
    output: "orthosnake/{region}/faa/{genome}.fasta"
    conda: "envs/scripts.yaml"
    shell:
        "name=$(basename {input});"
        "python scripts/GBfaa.py -gb  {input} > {output}"

rule make_ffn:
    input: "orthosnake/{region}/prokka/{genome}/{genome}.gbk"
    output: "orthosnake/{region}/ffn/{genome}.fasta"
    conda: "envs/scripts.yaml"
    shell:
        "name=$(basename {input});"
        "python scripts/GBffn.py -gb  {input} > {output}"


rule orthofinder:
    input: lambda wildcards:
        ["orthosnake/{0}/faa/{1}.fasta".format(wildcards.region, genome) for genome in get_selected_genomes(wildcards)]
    output: "orthosnake/{region}/Results/Orthogroups.txt"
    threads: 50
    params: folder="orthosnake/{region}"
    conda: "envs/ortho.yaml"
    shell:"bash scripts/run_orthofinder.sh {params.folder} {threads}"

checkpoint makeCoreOGfastasffn: 
    input:
        og="orthosnake/{region}/Results/Orthogroups.txt",
        ffns=lambda wildcards:
            ["orthosnake/{0}/ffn/{1}.fasta".format(wildcards.region, genome) for genome in get_selected_genomes(wildcards)]
    output: directory("orthosnake/{region}/Results/ortho/coreogs_nuc")
    params: folder="orthosnake/{region}"
    shell: 
        """
            set -eo pipefail
            folder="{params.folder}"
            mkdir -p $folder/tmp 
            cat $folder/ffn/*.fasta > $folder/tmp/all_genes_nuc.fasta
            mkdir -p $folder/Results/ortho/
            mkdir -p $folder/Results/ortho/coreogs_nuc
            perl scripts/splitToOg.pl $folder/Results/Orthogroups.txt $folder/tmp/all_genes_nuc.fasta $folder/Results/ortho/coreogs_nuc $folder/Results/Orthogroups_SingleCopyOrthologues.txt >> $folder/plog
        """ 

checkpoint makeCoreOGfastasfaa: 
    input:
        og="orthosnake/{region}/Results/Orthogroups.txt",
        faas=lambda wildcards:
            ["orthosnake/{0}/faa/{1}.fasta".format(wildcards.region, genome) for genome in get_selected_genomes(wildcards)]
    output: directory("orthosnake/{region}/Results/ortho/coreogs_prot")
    params: folder="orthosnake/{region}"
    shell: 
        """
            set -eo pipefail
            folder={params.folder}
            mkdir -p $folder/tmp 
            cat $folder/faa/*.fasta > $folder/tmp/all_genes_prot.fasta
            mkdir -p $folder/Results/ortho/
            mkdir -p $folder/Results/ortho/coreogs_prot
            perl scripts/splitToOg.pl $folder/Results/Orthogroups.txt $folder/tmp/all_genes_prot.fasta $folder/Results/ortho/coreogs_prot $folder/Results/Orthogroups_SingleCopyOrthologues.txt >> $folder/plog
        """ 

rule align_core_prot:
    input:
        "orthosnake/{region}/Results/ortho/coreogs_prot/{og}.fasta"
    output:
        "orthosnake/{region}/Results/ortho/coreogs_aligned_prot/{og}.fasta"
    shell:
        "helpers/./muscle -in {input} -out {output} -quiet"

rule pal2nal:
    input:
        prot="orthosnake/{region}/Results/ortho/coreogs_aligned_prot/{og}.fasta",
        nuc="orthosnake/{region}/Results/ortho/coreogs_nuc/{og}.fasta"
    output: "orthosnake/{region}/Results/ortho/coreogs_aligned_nuc/{og}.fasta"
    shell:
        "perl helpers/pal2nal.pl {input.prot} {input.nuc} -output fasta > {output}"

def aggregate_core_og(wildcards):
    checkpoint_output = checkpoints.makeCoreOGfastasffn.get(**wildcards).output[0]
    ogs = expand("orthosnake/{region}/Results/ortho/coreogs_aligned_nuc/{og}.fasta",region=wildcards.region,
           og=glob_wildcards(os.path.join(checkpoint_output, "{og}.fasta")).og)
    return ogs

rule cat_core:
    input: aggregate_core_og
    output: "orthosnake/{region}/tmp/coreogaligned.fasta" 
    conda: "envs/scripts_tree.yaml"
    params: folder="orthosnake/{region}"
    shell:
        "perl scripts/concatenate_core.pl {params.folder}/Results/ortho/coreogs_aligned_nuc {output}" 

rule tree_for_core:
    input: "orthosnake/{region}/tmp/coreogaligned.fasta" 
    threads: 20
    params: folder="orthosnake/{region}"
    output: "orthosnake/{region}/Results/coreogs_nucleotide.treefile"
    shell:
        "iqtree -s {input} --seqtype CODON -T AUTO --threads-max {threads} --prefix {params.folder}/Results/coreogs_nucleotide -redo -m MFP"



#########################
## TREE OF REGION OG

rule og_in_region:
    input: "orthosnake/{region}/Results/Orthogroups.txt"
    output: "tmp/og_in_region/{region}/og_list"
    conda: "envs/env.yaml"
    shell: "Rscript helpers/genes_region_overlap.r {wildcards.region}"
import os
import re
import glob
import random

results = []

genomes,=glob_wildcards("fna/{genome}.fna")
#genomes = genomes[0:60]

regions,=glob_wildcards("regions/{region}.fasta")
results.append(expand("tree_region/{region}/{region}.treefile", region =  regions))

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
    params: min_coverage=90, max_genomes = 100
    run:
        region = wildcards.region
        os.mkdir(output.outdir)
        genomes_with_regions_list = []
        for cov_file in input:
            if(os.path.getsize(cov_file) < 1): continue
            with open(cov_file, "r") as infile:
                line = infile.readline()
                print(line)
                [genome, coverage] = re.split(r'\t+', line)
                if(int(coverage) > params.min_coverage): 
                    genomes_with_regions_list.append(genome)
        if(len(genomes_with_regions_list) > params.max_genomes):
            genomes_with_regions_list = random.sample(genomes_with_regions_list, params.max_genomes)                        
        for genome in genomes_with_regions_list:
            os.symlink("../../fna/"+genome+".fna", "genomes_with_regions/"+region+"/"+genome+".fna")

rule nucmer_region_to_genomes:
    input: genome="genomes_with_regions/{region}/{genome}.fna",
            region="regions/{region}.fasta"
    output: "tmp/nucmer/{region}/{genome}.coord"
    params: prefix="tmp/nucmer/{region}/{genome}"
    conda: "envs/env.yaml"
    threads: 1
    shell:
        """
            nucmer {input.genome} {input.region} -p {params.prefix} --mum -L 1000 -t {threads}
            show-coords -HlcT {params.prefix}.delta > {output}
        """

def aggregate_by_selected_genomes(wildcards):
    selected_genomes_folder = checkpoints.select_genomes_with_regions.get(**wildcards).output[0]
    result = expand("tmp/nucmer/{region}/{genome}.coord",
           region=wildcards.region,
           genome=glob_wildcards(os.path.join( selected_genomes_folder, "{genome}.fna")).genome)
    return result

rule mauve_of_region: 
    input: aggregate_by_selected_genomes
    output: "tmp/result/{region}/out"
    shell:
        """ 
        cat {input} > {output}
        """

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
    output: "orthosnake/{region}/Results/Orthogroups.txt", "orthosnake/{region}/Results/Orthogroups_SingleCopyOrthologues.txt"
    threads: 50
    params: folder="orthosnake/{region}"
    conda: "envs/ortho.yaml"
    shell:"bash scripts/run_orthofinder.sh {params.folder} {threads}"

rule cat_genes:
    input: ffns=lambda wildcards:
            ["orthosnake/{0}/ffn/{1}.fasta".format(wildcards.region, genome) for genome in get_selected_genomes(wildcards)]
    output: "orthosnake/{region}/tmp/all_genes_nuc.fasta"
    shell:
        """
        mkdir -p orthosnake/{wildcards.region}/tmp 
        cat orthosnake/{wildcards.region}/ffn/*.fasta > {output}
        """

checkpoint makeCoreOGfastasffn: 
    input:
        og="orthosnake/{region}/Results/Orthogroups.txt",
        all_genes="orthosnake/{region}/tmp/all_genes_nuc.fasta",
        core_og="orthosnake/{region}/Results/Orthogroups_SingleCopyOrthologues.txt"
    output: directory("orthosnake/{region}/Results/coreogs_nuc")
    shell: 
        """
            python scripts/splitToOg.py {input.og} {input.all_genes} {input.core_og} {output}
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
    ## add iqtree to conda
    output: "orthosnake/{region}/Results/coreogs_nucleotide.treefile"
    shell:
        "iqtree -s {input} --seqtype CODON -T AUTO --threads-max {threads} --prefix {params.folder}/Results/coreogs_nucleotide -redo -m MFP"

#============================================ REGION TREE --

rule og_in_region:
    input:  lambda wildcards:
        ["blast/{0}/{1}.blast".format(wildcards.region, get_selected_genomes(wildcards)[0]),
        "orthosnake/{0}/faa/{1}.fasta".format(wildcards.region, get_selected_genomes(wildcards)[0])]
    output: "tmp/og_in_region/{region}/og_list"
    shell:
        "Rscript helpers/genes_region_overlap.r {wildcards.region} {output}"

rule og_in_region_sequence:
    input: oglist="tmp/og_in_region/{region}/og_list", ogaligned="orthosnake/{region}/tmp/coreogaligned.fasta"
    output: "tmp/og_in_region/{region}/{region}.fasta"
    params: folder="orthosnake/{region}"
    shell: 
        """
            while read f
            do 
                OG_FASTA="orthosnake/{wildcards.region}/Results/ortho/coreogs_aligned_nuc/$f.fasta"
                if [ -f $OG_FASTA ]; then  
                    cp $OG_FASTA tmp/og_in_region/{wildcards.region}/$f.fasta
                fi    
            done < {input.oglist}
            perl scripts/concatenate_core.pl tmp/og_in_region/{wildcards.region} {output}
        """

rule tree_for_region:
    input: "tmp/og_in_region/{region}/{region}.fasta"
    output: "tree_region/{region}/{region}.treefile"
    params: folder="tree_region/{region}/"
    ## add iqtree to conda
    threads: 20
    shell:
        "iqtree -s {input} --seqtype CODON -T AUTO --threads-max {threads} --prefix tree_region/{wildcards.region}/{wildcards.region} -redo -m MFP"

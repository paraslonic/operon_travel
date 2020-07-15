import os
import re
import glob

results = []

genomes,=glob_wildcards("fna/{genome}.fna")
#genomes = genomes[0:60]

#results.append(expand("coverage/pdu/{genome}.cov", genome=genomes))
#results.append(expand("blast/{genome}.blast", genome=genomes))

#results.append("genomes_with_regions/pdu")
results.append("tmp/pdu/mauve.backbone")

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

rule find_region_in_selected_genome:
    input:  blast="blast/{region}/{genome}.blast",
            genome="genomes_with_regions/{region}/{genome}.fna"
    output: bed="tmp/subject_regions/{region}/{genome}.bed",
            fasta="tmp/subject_regions/{region}/{genome}.fasta",
    conda: "envs/env.yaml"
    shell: 
        """
            Rscript helpers/regionFromBlast.r {input.blast} {output.bed}
            bedtools getfasta -fi {input.genome} -fo {output.fasta} -bed {output.bed}
            #touch {output.bed}
        """

def aggregate_by_selected_genomes(wildcards):
    selected_genomes_folder = checkpoints.select_genomes_with_regions.get(**wildcards).output[0]
    result = expand("tmp/subject_regions/{region}/{genome}.fasta",
           region=wildcards.region,
           genome=glob_wildcards(os.path.join( selected_genomes_folder, "{genome}.fna")).genome)
    return result

rule get_list:
    input: subject_regions=aggregate_by_selected_genomes,
            query_region="regions/{region}.fasta"
    output: backbone="tmp/{region}/mauve.backbone",
            xmfa="tmp/{region}/mauve.xmfa"
    shell:
        """ 
        #echo {input} > {output}
        progressiveMauve {input.query_region} {input.subject_regions} --output {output.xmfa} --backbone-output {output.backbone}
        """
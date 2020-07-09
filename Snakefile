import os
import re

results = []

genomes,=glob_wildcards("fna/{genome}.fna")
genomes = genomes[0:120]

results.append(expand("coverage/pdu/{genome}.cov", genome=genomes))
results.append("genomes_with_regions/pdu")

print("Genomes: ", len(genomes))

rule all:
    input: results

rule blast:
    input: region="regions/{region}.fasta",genome="fna/{genome}.fna"
    output: "blast/{region}/{genome}.blast"
    shell: 
        """
        blastn -query {input.region} -subject {input.genome} -outfmt '6 std qcovs' -out {output}
        """

rule calculate_coverage:
    input: "blast/{region}/{genome}.blast"
    output: "coverage/{region}/{genome}.cov"
    shell: 
        """
        printf "{wildcards.genome}\t" > {output}
        if [ -s "{input}" ] 
        then
            head -1  {input} | awk '{{print $NF}}'  >> {output}
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


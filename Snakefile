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
	shell:	"python scripts/cropHeader.py -input  {input} -out {output} -n 20"
		
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
	input:	"orthosnake/{region}/prokka/{genome}/{genome}.gbk"
	output: "orthosnake/{region}/faa/{genome}.fasta"
	conda: "envs/scripts.yaml"
	shell:
		"name=$(basename {input});"
		"python scripts/GBfaa.py -gb  {input} > {output}"

rule make_ffn:
	input:	"orthosnake/{region}/prokka/{genome}/{genome}.gbk"
	output: "orthosnake/{region}/ffn/{genome}.fasta"
	conda: "envs/scripts.yaml"
	shell:
		"name=$(basename {input});"
		"python scripts/GBffn.py -gb  {input} > {output}"

rule orthofinder:
	input: 
		expand("faa/{genome}.fasta", genome = get_selected_genomes(wildcards))
	output:"orthosnake/{region}/Results/Orthogroups.txt"
	threads: 50
    params: dir="orthosnake/{region}/faa"
	conda: "envs/ortho.yaml"
	log: "log_of.txt"
	shell:
		"bash scripts/run_orthofinder.sh {params.dir} {threads} > {log}" #### MODIFY SET INPUT FOLDER

checkpoint makeCoreOGfastas: 
	input:
		og="orthosnake/{region}/Results/Orthogroups.txt",
		ffns=expand("orthosnake/{region}/ffn/{qu}.fasta", qu=GENOMES)
	output:
		coreog=directory("orthosnake/{region}/Results/ortho/coreogs_prot")
	shell:
		"""
		mkdir -p tmp #### MODIFY !!!
		cat ffn/*.fasta > tmp/all_genes_nuc.fasta
		cat faa/*.fasta > tmp/all_genes_prot.fasta
        mkdir -p Results/ortho/
		mkdir -p Results/ortho/coreogs_nuc Results/ortho/coreogs_prot
		perl scripts/splitToOg.pl Results/Orthogroups.txt tmp/all_genes_nuc.fasta Results/ortho/coreogs_nuc Results/Orthogroups_SingleCopyOrthologues.txt
		perl scripts/splitToOg.pl Results/Orthogroups.txt tmp/all_genes_prot.fasta Results/ortho/coreogs_prot Results/Orthogroups_SingleCopyOrthologues.txt
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
	output:
		"orthosnake/{region}/Results/ortho/coreogs_aligned_nuc/{og}.fasta"
	shell:
		"perl helpers/pal2nal.pl {input.prot} {input.nuc} -output fasta > {output}"

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.makeCoreOGfastas.get(**wildcards).output[0]
    ogs = expand("orthosnake/{region}/Results/ortho/coreogs_aligned_nuc/{og}.fasta",
           og=glob_wildcards(os.path.join(checkpoint_output, "{og}.fasta")).og)
    return ogs

rule cat_core:
    input: aggregate_input
    output: "orthosnake/{region}/tmp/coreogaligned.fasta" 
    conda: "envs/scripts_tree.yaml"
    shell:
        "echo {input};"
        "perl scripts/concatenate_core.pl Results/ortho/coreogs_aligned_nuc {output}" ### MODIFY

rule tree_for_core:
	input: rules.cat_core.output
	threads: 20
	output:
		"orthosnake/{region}/Results/coreogs_nucleotide.treefile"
	shell:
		"iqtree -s {input} -nt {threads} -pre Results/coreogs_nucleotide -redo -m MFP"
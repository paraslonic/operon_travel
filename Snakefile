results = []

genomes,=glob_wildcards("fna/{genome}.fna")
genomes = genomes[0:20]

results.append(expand("coverage/pdu/{genome}.cov", genome=genomes))
results.append("genomes_with_regions/pdu/GCF_012934785.1_ASM1293478v1_genomic.fna")

print("Genomes: ", len(genomes))

rule all:
    input: results

checkpoint select_genomes_with_regions:
    input: coverage="coverage/{region}/{genome}.cov"
    output: "genomes_with_regions/{region}/{genome}.fna"
    shell: 
        """
        awk '{{ if($2>90) {{ system("ln -s $(pwd)/fna/{wildcards.genome}.fna {output}") }} }}' {input}
        #touch {output}
        """

rule calculate_coverage:
    input: "blast/{region}/{genome}.blast"
    output: "coverage/{region}/{genome}.cov"
    shell: 
        """
        printf "{wildcards.genome}\t" > {output}
        head -1 {input} | awk '{{print $NF}}'  >> {output}
        """

rule blast:
    input: region="regions/{region}.fasta",genome="fna/{genome}.fna"
    output: "blast/{region}/{genome}.blast"
    shell: 
        """
        blastn -query {input.region} -subject {input.genome} -outfmt '6 std qcovs' -out {output}
        """
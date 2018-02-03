
rule all:
    input: roary

rule esearch:
    input: 'batch_entrez.txt'
    output: 'ncbi/{sample}.fasta'
    shell: 'while read -r line; do'
    esearch :

rule gbk_to_gff3:

    input: '{sample}.gbk'
    output: 'gff3/{sample}.gff'
    shell: 'bp_genbank2gff3.pl -o'

rule roary:
    input: 'gff3/{sample}.gff'
    output: 'roary/'
    shell: 'roary {input}'


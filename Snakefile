READS = ["1", "2"]

rule download_data:
    output:
        ref="data/GCA_007002865.1_ASM700286v1_genomic.fna.gz",
        r1="data/SRR12827773_1.fastq",
        r2="data/SRR12827773_2.fastq"
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Agrobacterium_tumefaciens/latest_assembly_versions/GCA_007002865.1_ASM700286v1/GCA_007002865.1_ASM700286v1_genomic.fna.gz -P data 
        fastq-dump SRR12827773 --split-files -O data
        """

rule fastqc:
    input:
        r1=rules.download_data.output.r1,
        r2=rules.download_data.output.r2
    output:
        expand("fastqc/SRR12827773_{r}_fastqc.html", r=READS),
        expand("fastqc/SRR12827773_{r}_fastqc/summary.txt", r=READS)
    shell:
        """
        fastqc -o fastqc {input.r1} {input.r2}
        cd fastqc
        unzip SRR12827773_1_fastqc.zip
        unzip SRR12827773_2_fastqc.zip
        rm -f SRR12827773_*_fastqc.zip
        cd ../
        """

rule spades:
    input:
        r1=rules.download_data.output.r1,
        r2=rules.download_data.output.r2
    params:
        outdir="spades"
    output:
        "{outdir}/scaffolds.fasta",
        "{outdir}/contigs.fasta"
    shell:
        "spades.py -o spades -1 {input.r1} -2 {input.r2}"
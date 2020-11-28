rule all:
    input:
        "PAO1.AiiA-lactonase.variant_effects.png",
        "PAO1.AiiA-lactonase.qs_variants.tsv",
        "PAO1.AiiA-lactonase.qs_variants.png"

rule download_ref:
    output:
        ref="data/GCF_000006765.1_ASM676v1_genomic.fna"
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1_genomic.fna.gz -P data
        gunzip data/GCF_000006765.1_ASM676v1_genomic.fna.gz
        """

rule download_reads:
    output:
        r1="data/SRR12820667_1.fastq",
        r2="data/SRR12820667_2.fastq"
    shell:
        """
        fastq-dump SRR12820667 --split-files -O data
        """

rule bwa_mem:
    input:
        ref=rules.download_ref.output.ref,
        r1=rules.download_reads.output.r1,
        r2=rules.download_reads.output.r2
    output:
        bam="alignments.sorted.bam"
    threads: 8
    shell:
        """
        bwa index {input.ref}
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -u -F 4 -q 30 -@ {threads} | samtools sort -O BAM  -o {output.bam} -@ {threads}
        """

rule samtools_faidx:
    input:
        ref=rules.download_ref.output.ref
    output:
        idx="data/GCF_000006765.1_ASM676v1_genomic.fna.fai"
    shell:
        """
        samtools faidx {input.ref}
        """

rule call_variants:
    input:
        bam=rules.bwa_mem.output.bam,
        ref=rules.download_ref.output.ref,
        idx=rules.samtools_faidx.output.idx
    output:
        vcf="variants.named.vcf.gz"
    threads: 8
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} {input.bam} | bcftools call -Ou -mv | bcftools filter -e 'QUAL<30 | DP<20' --threads {threads} -Oz > variants.vcf.gz
        gunzip -c variants.vcf.gz | sed -E 's/^NC_002516.2/Chromosome/g' | gzip > {output.vcf}
        """

rule annotate_variants:
    input:
        vcf=rules.call_variants.output.vcf
    output:
        snpEff="variants.named.snpEff.vcf",
        genes="variants.snpEff.summary.genes.txt",
        csv="variants.snpEff.summary.csv"
    shell:
        """
        snpEff ann -csvStats {output.csv} -geneId Pseudomonas_aeruginosa_pao1 {input.vcf} > {output.snpEff}
        """

rule parse_variants:
    input:
        genes=rules.annotate_variants.output.genes,
        ref_qs="PAO1_quorum_sensing_genes.tsv"
    output:
        var_plot="PAO1.AiiA-lactonase.variant_effects.png",
        qs_tsv="PAO1.AiiA-lactonase.qs_variants.tsv",
        qs_plot="PAO1.AiiA-lactonase.qs_variants.png"
    shell:
        """
        ./parse_variants.py {input.genes} {input.ref_qs}
        """

rule dag:
    output:
        dag="snakemake_workflow_dag.png"
    shell:
        """
        snakemake parse_variants --dag | dot -Tpng > {output.dag}
        """
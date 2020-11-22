READS = ["1", "2"] 
# get # threads as a param from the command line

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
    shell:
        """
        bwa index {input.ref}
        bwa mem {input.ref} {input.r1} {input.r2} | samtools view -u -F 4 -q 30 | samtools sort -O BAM  -o {output.bam}
        """

rule bwa_index:
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
        idx=rules.bwa_index.output.idx
    output:
        vcf="variants.named.vcf.gz"
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} {input.bam} | bcftools call -Ou -mv | bcftools filter -e 'QUAL<30 | DP<20' -Ou | sed -r s/NC_002516.2/Chromosome/g | gzip > variants.named.vcf.gz
        """

rule annotate_variants:
    input:
        vcf=rules.call_variants.output.vcf
    output:
        snpEff="variants.named.snpEff.vcf",
        csv="variants.snpEff.summary.csv"
    shell:
        """
        snpEff ann -csvStats {output.csv} Pseudomonas_aeruginosa_pao1 {input.vcf} > {output.snpEff}
        """

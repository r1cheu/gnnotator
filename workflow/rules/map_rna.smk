rule hisat2_index:
    input:
        masked_fa="results/mask_repeat/{species}.masked.fa",
    output:
        index_file=multiext(
            "results/map_rna/{species}/{species}",
            ".1.ht2",
            ".2.ht2",
            ".3.ht2",
            ".4.ht2",
            ".5.ht2",
            ".6.ht2",
            ".7.ht2",
            ".8.ht2",
        ),
    log:
        out="logs/map_rna/{species}/hisat2_index.out",
        err="logs/map_rna/{species}/hisat2_index.err",
    params:
        prefix="results/map_rna/{species}/{species}",
    conda:
        "../envs/map_rna.yml"
    shell:
        """
        hisat2-build {input.masked_fa} {params.prefix} 2> {log.err} 1> {log.out}
        """


rule hisat2_align:
    input:
        index=rules.hisat2_index.output.index_file,
        r1="data/rnaseq/{species}_{tissue}_R1.fastq",
        r2="data/rnaseq/{species}_{tissue}_R2.fastq",
    output:
        sorted_bam="results/map_rna/{species}/{tissue}_sorted.bam",
    params:
        prefix="results/map_rna/{species}/{tissue}",
        index_prefix="results/map_rna/{species}/{species}",
    log:
        out="logs/map_rna/{species}/{tissue}_hisat2.out",
        err="logs/map_rna/{species}/{tissue}_hisat2.err",
    conda:
        "../envs/map_rna.yml"
    shell:
        """
        hisat2 -x {params.index_prefix} -1 {input.r1} -2 {input.r2} --dta --rg-id {wildcards.tissue} --rg "SM:{wildcards.tissue},LB:lib1,PL:ILLUMINA" -S {params.prefix}.sam -p 8 2> {log.err} 1> {log.out}

        samtools view -@ 8 -bS {params.prefix}.sam | samtools sort -@ 8 -o {output.sorted_bam} - 2 >{log.err} 1> {log.out}

        rm {params.prefix}.sam
        """


rule merge_bams:
    input:
        bams=expand(
            "results/map_rna/{species}/{tissue}_sorted.bam",
            species="{species}",
            tissue=["fringe", "leaf", "root", "seedling"],
        ),
    output:
        merged_bam="results/map_rna/{species}/Merge.bam",
    conda:
        "../envs/map_rna.yml"
    log:
        err="logs/map_rna/{species}/merge.err",
    shell:
        """
        samtools merge -@ 8 -o {output.merged_bam} {input.bams} 2> {log.err}
        """


rule stringtie:
    input:
        bam="results/map_rna/{species}/Merge.bam",
    output:
        merged_gtf="results/map_rna/{species}/Merge.gtf",
    conda:
        "../envs/map_rna.yml"
    shell:
        """
        stringtie -o {output.merged_gtf} -p 8 {input.bam}
        """

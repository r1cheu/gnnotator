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
        "logs/map_rna/{species}/hisat2_index.log",
    params:
        prefix=lambda w, output: output.index_file[0].replace(".1.ht2", ""),
    conda:
        "../envs/map_rna.yml"
    shell:
        """
        hisat2-build {input.masked_fa} {params.prefix} &>{log}
        """


rule hisat2_align:
    input:
        index=rules.hisat2_index.output.index_file,
        r1="data/rnaseq/{species}_{tissue}_R1.fastq",
        r2="data/rnaseq/{species}_{tissue}_R2.fastq",
    output:
        sorted_bam="results/map_rna/{species}/{tissue}_sorted.bam",
    params:
        prefix=lambda w, output: output.sorted_bam.replace("_sorted.bam", ""),
        index_prefix=lambda w, input: input.index[0].replace(".1.ht2", ""),
    log:
        "logs/map_rna/{species}/{tissue}_hisat2.log",
    conda:
        "../envs/map_rna.yml"
    threads: 8
    resources:
        cpus_per_task=threads,
    shell:
        """
        hisat2 -x {params.index_prefix} -1 {input.r1} -2 {input.r2} --dta --rg-id {wildcards.tissue} --rg "SM:{wildcards.tissue},LB:lib1,PL:ILLUMINA" -S {params.prefix}.sam -p {threads} &>{log}
        samtools view -@ {threads} -bS {params.prefix}.sam | samtools sort -@ {threads} -o {output.sorted_bam} - &>>{log}
        rm {params.prefix}.sam
        """


rule merge_bams:
    input:
        bams=expand(
            "results/map_rna/{species}/{tissue}_sorted.bam",
            species="{species}",
            tissue=config["rna_seq_tissues"],
        ),
    output:
        merged_bam="results/map_rna/{species}/Merge.bam",
    conda:
        "../envs/map_rna.yml"
    log:
        "logs/map_rna/{species}/merge.log",
    threads: 8
    resources:
        cpus_per_task=threads,
    shell:
        """
        samtools merge -@ {threads} -o {output.merged_bam} {input.bams} &>{log}
        """


rule stringtie:
    input:
        bam="results/map_rna/{species}/Merge.bam",
    output:
        merged_gtf="results/map_rna/{species}/Merge.gtf",
    conda:
        "../envs/map_rna.yml"
    log:
        "results/map_rna/{species}/stringtie.log",
    threads: 8
    resources:
        cpus_per_task=threads,
    shell:
        """
        stringtie -o {output.merged_gtf} -p {threads} {input.bam} &>{log}
        """

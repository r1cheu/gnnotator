rule trinity_genome_guided:
    input:
        bam="results/map_rna/{species}/Merge.bam",
    output:
        fa="results/trinity/{species}/trinity_genome_guided.fasta",
    log:
        "logs/trinity/{species}/genome_guided.log",
    params:
        dir=lambda w, output: output.fa.replace(".fasta", ""),
    conda:
        "../envs/trinity.yml"
    threads: 32
    resources:
        mem_mb=102400,
        cpus_per_task=threads,
    shell:
        """
        Trinity --genome_guided_bam {input.bam} \
            --genome_guided_max_intron 10000 \
            --max_memory 200G \
            --CPU {threads} \
            --output {params.dir} &>{log}
        mv {params.dir}/Trinity-GG.fasta {output.fa}
        rm -rf {params.dir}
        """


rule trinity_de_novo:
    input:
        r1=expand(
            "data/rnaseq/{species}_{tissue}_R1.fastq",
            species="{species}",
            tissue=["fringe", "leaf", "root", "seedling"],
        ),
        r2=expand(
            "data/rnaseq/{species}_{tissue}_R2.fastq",
            species="{species}",
            tissue=["fringe", "leaf", "root", "seedling"],
        ),
    output:
        fa="results/trinity/{species}/trinity_de_novo.Trinity.fasta",
    log:
        "logs/trinity/{species}/denovo.log",
    conda:
        "../envs/trinity.yml"
    params:
        dir=lambda w, output: output.fa.replace(".Trinity.fasta", ""),
        left_reads=lambda wildcards, input: ",".join(input.r1),
        right_reads=lambda wildcards, input: ",".join(input.r2),
    threads: 32
    resources:
        mem_mb=102400,
        cpus_per_task=threads,
    shell:
        """
        Trinity --seqType fq \
            --max_memory 100G \
            --left {params.left_reads} \
            --right {params.right_reads} \
            --CPU {threads} \
            --output {params.dir} &> {log}

        rm -rf {params.dir}
        """


rule combine_trinity_outputs:
    input:
        denovo=rules.trinity_de_novo.output.fa,
        genome_guided=rules.trinity_genome_guided.output.fa,
    output:
        combined="results/trinity/{species}/{species}.fa",
    log:
        "logs/trinity/{species}/combine.log",
    conda:
        "../envs/base.yml"
    shell:
        """
        cat {input.denovo} {input.genome_guided} > {output.combined}
        echo "Combined Trinity outputs into {output.combined}" &> {log}
        """

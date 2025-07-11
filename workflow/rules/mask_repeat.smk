rule mask_repeat:
    input:
        genome="data/assembly/{species}.fa",
    output:
        masked="results/mask_repeat/{species}.masked.fa",
        tbl="results/mask_repeat/{species}.masked.fa.tbl",
    container:
        container_image("leizzzz/repeatmasker:4.0.6")
    message:
        "Masking repeats in {wildcards.species} assembly"
    log:
        "logs/mask_repeat/{species}.log",
    threads: 8
    resources:
        cpus_per_task=threads,
    shell:
        """
        RepeatMasker -engine rmblast -pa {threads} -nolow -species rice -dir results/mask_repeat/ {input.genome} &>{log}

        mv results/mask_repeat/{wildcards.species}.fa.masked {output.masked}

        mv results/mask_repeat/{wildcards.species}.fa.tbl {output.tbl} 
        """

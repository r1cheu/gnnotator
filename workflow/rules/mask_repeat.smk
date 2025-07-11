import os
from pathlib import Path


rule mask_repeat:
    input:
        genome="data/assembly/{species}.fa",
    output:
        masked="results/mask_repeat/{species}/{species}.masked.fa",
        tbl="results/mask_repeat/{species}/{species}.masked.fa.tbl",
    container:
        container_image("leizzzz/repeatmasker:4.0.6")
    message:
        "Masking repeats in {wildcards.species} assembly"
    log:
        "logs/mask_repeat/{species}.log",
    params:
        genome=lambda w, input: Path(input.genome).resolve(),
        dir=lambda w, output: Path(output.masked).parent.resolve(),
        masked=lambda w, output: Path(output.masked).resolve(),
        tbl=lambda w, output: Path(output.tbl).resolve(),
    threads: 8
    resources:
        cpus_per_task=threads,
    shell:
        """
        cd {params.dir}
        micromamba run -n base RepeatMasker -engine rmblast -pa {threads} -nolow -species rice -dir {params.dir} {params.genome} &>{log}

        mv {wildcards.species}.fa.masked {params.masked}
        mv {wildcards.species}.fa.tbl {params.tbl} 
        """

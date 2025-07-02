rule install_repeatmasker_library:
    input:
        libraries="resources/repeatmaskerlibraries-20150807.tar.gz",
    output:
        touch("results/repeatmaskerlibraries-20150807.done"),
    conda:
        "../envs/mask_repeat.yml"
    log:
        "logs/install_repeatmasker_library.log",
    shell:
        """
        RM_EXECUTABLE=$(which RepeatMasker)
        TARGET_DIR=$(dirname "${{RM_EXECUTABLE}}")/../share/RepeatMasker/

        echo "Found RepeatMasker executable at: ${{RM_EXECUTABLE}}" &>{log}
        echo "Target directory for RepeatMasker libraries: ${{TARGET_DIR}}" &>>{log}

        echo "Extracting RepeatMasker libraries ..." &>>{log}
        tar -xzf {input.libraries} -C ${{TARGET_DIR}}
        echo "RepeatMasker libraries installed successfully." &>>{log}
        """


rule mask_repeat:
    input:
        libraries_installed="results/repeatmaskerlibraries-20150807.done",
        genome="data/assembly/{species}.fa",
    output:
        masked="results/mask_repeat/{species}.masked.fa",
        tbl="results/mask_repeat/{species}.masked.fa.tbl",
    conda:
        "../envs/mask_repeat.yml"
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

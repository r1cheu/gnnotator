
rule install_repeatmasker_library:
    input:
        libraries="resource/repeatmaskerlibraries-20150807.tar.gz",
    output:
        touch("results/repeatmaskerlibraries-20150807.done"),
    conda:
        "../envs/mask_repeat.yml"
    shell:
        """
        RM_EXECUTABKE=$(which RepeatMasker)
        TARGET_DIR=$(dirname "${{RM_EXECUTABLE}}")/../share/RepeatMasker/Library

        echo "Found RepeatMasker executable at: ${{RM_EXECUTABLE}}"
        echo "Target directory for RepeatMasker libraries: ${{TARGET_DIR}}"

        mkdir -p ${{TARGET_DIR}}

        echo "Extracting RepeatMasker libraries to ${{TARGET_DIR}}"

        tar -xzf {input.libraries} -C ${{TARGET_DIR}}
        """


rule mask_repeat:
    input:
        libraries_installed="results/repeatmaskerlibraries-20150807.done",
        genome="data/{spcies}.fa",
    output:
        masked="results/mask_repeat/{spcies}.masked.fa",
        tbl="results/mask_repeat/{spcies}.masked.fa.tbl",
    conda:
        "../envs/mask_repeat.yml"
    message:
        "Masking repeats in {wildcards.spcies} assembly"
    log:
        "results/mask_repeat/{spcies}.log",
    shell:
        """
        /usr/bin/bash
        RepeatMasker -engine rmblast -pa 20 -nolow -species rice -dir results/mask_repeat/{wildcards.spcies} {input.genome} > {log} 2>&1 \
        mv results/mask_repeat/{wildcards.spcies}/{wildcards.spcies}.fa.masked {output.masked} \
        mv results/mask_repeat/{wildcards.spcies}/{wildcards.spcies}.fa.tbl {output.tbl} 
        """

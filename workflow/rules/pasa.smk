from pathlib import Path


rule create_pasa_config:
    output:
        pasa_config="results/pasa/{species}/pasa.alignAssembly.txt",
    params:
        conf=config["pasa_config"],
    conda:
        "../envs/base.yml"
    log:
        "logs/pasa/{species}/create_pasa_config.log",
    script:
        "../scripts/create_pasa_config.py"


rule pasa:
    input:
        pasa_config="results/pasa/{species}/pasa.alignAssembly.txt",
        transcripts="results/trinity/{species}/{species}.fa",
        masked_genome="results/mask_repeat/{species}.masked.fa",
        trans_gtf="results/map_rna/{species}/Merge.gtf",
    output:
        assemblies_gff="results/pasa/{species}/{species}_DB_pasa.sqlite.pasa_assemblies.gff3",
        blat_gff="results/pasa/{species}/{species}_DB_pasa.sqlite.valid_blat_alignments.gff3",
        gmap_gff="results/pasa/{species}/{species}_DB_pasa.sqlite.valid_gmap_alignments.gff3",
        assemblies="results/pasa/{species}/{species}_DB_pasa.sqlite.assemblies.fasta",
    params:
        pasa_dir=lambda w, input: os.path.dirname(input.pasa_config),
        abs_pasa_config=lambda w, input: Path(input.pasa_config).resolve(),
        abs_masked_genome=lambda w, input: Path(input.masked_genome).resolve(),
        abs_transcripts=lambda w, input: Path(input.transcripts).resolve(),
        abs_trans_gtf=lambda w, input: Path(input.trans_gtf).resolve(),
        abs_log=lambda w: Path(f"logs/pasa/{w.species}/pasa.log").resolve(),
    container:
        container_image("pasapipeline/pasapipeline:2.5.3")
    log:
        "logs/pasa/{species}/pasa.log",
    threads: 16
    resources:
        mem_mb=10240,
        cpus_per_task=threads,
    shell:
        """
        cd {params.pasa_dir}
        /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
            -c {params.abs_pasa_config} \
            -C -R \
            -g {params.abs_masked_genome} \
            -t {params.abs_transcripts} \
            --trans_gtf {params.abs_trans_gtf} \
            --ALIGNERS blat,gmap \
            --CPU {threads} &>{params.abs_log}
        """


rule pasa_dbi:
    input:
        fa=rules.pasa.output.assemblies,
        gff=rules.pasa.output.assemblies_gff,
    output:
        pasa_dbi="results/pasa/{species}/{species}_DB_pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3",
    params:
        pasa_dir=lambda w, output: os.path.dirname(output.pasa_dbi),
        abs_fa=lambda wildcards, input: Path(input.fa).resolve(),
        abs_gff=lambda wildcards, input: Path(input.gff).resolve(),
        abs_log=lambda w: Path(f"logs/pasa/{w.species}/pasa_dbi.log").resolve(),
    log:
        "logs/pasa/{species}/pasa_dbi.log",
    container:
        container_image("pasapipeline/pasapipeline:2.5.3")
    resources:
        mem_mb=10240,
    shell:
        """
        cd {params.pasa_dir}
        /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta {params.abs_fa} --pasa_transcripts_gff3 {params.abs_gff} &> {params.abs_log}
        """


rule cat_all_gff:
    input:
        gff=rules.pasa.output.assemblies_gff,
        blat_gff=rules.pasa.output.blat_gff,
        gmap_gff=rules.pasa.output.gmap_gff,
        transdecoder_gff=rules.pasa_dbi.output.pasa_dbi,
    output:
        all_gff="results/pasa/{species}/{species}.gff3",
    conda:
        "../envs/base.yml"
    log:
        "logs/pasa/{species}/cat_all_gff.log",
    shell:
        """
        cat {input.gff} {input.blat_gff} {input.gmap_gff} {input.transdecoder_gff} > {output.all_gff}
        echo "Created combined GFF3 file with PASA assemblies, BLAT alignments, GMAP alignments, and TransDecoder annotations." &> {log}
        """

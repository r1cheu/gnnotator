from pathlib import Path


rule create_pasa_config:
    output:
        pasa_config="results/pasa/{species}/pasa.alignAssembly.txt",
    params:
        conf=config["pasa_config"],
    run:
        db_name = f"/tmp/{wildcards.species}_DB_pasa.sqlite"
        script_params = params.conf["params"]

        with open(output.pasa_config, "w") as f:
            f.write(f"DATABASE={db_name}\n")
            for script_name, params_dict in script_params.items():
                for param_key, param_val in params_dict.items():
                    f.write(f"{script_name}:{param_key}={param_val}\n")


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
        pasa_dir="results/pasa/{species}",
        abs_pasa_config=lambda wildcards, input: Path(input.pasa_config).resolve(),
        abs_masked_genome=lambda wildcards, input: Path(input.masked_genome).resolve(),
        abs_transcripts=lambda wildcards, input: Path(input.transcripts).resolve(),
        abs_trans_gtf=lambda wildcards, input: Path(input.trans_gtf).resolve(),
        log=lambda wildcards, input: Path(
            f"logs/pasa/{wildcards.species}/pasa.log"
        ).resolve(),
    container:
        "docker://pasapipeline/pasapipeline:2.5.3"
    log:
        "logs/pasa/{species}/pasa.log",
    threads: 32
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
            --CPU {threads} &> {params.log}
        """


rule pasa_dbi:
    input:
        fa=rules.pasa.output.assemblies,
        gff=rules.pasa.output.assemblies_gff,
    output:
        pasa_dbi="results/pasa/{species}/{species}_DB_pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3",
    params:
        pasa_dir="results/pasa/{species}",
        abs_fa=lambda wildcards, input: Path(input.fa).resolve(),
        abs_gff=lambda wildcards, input: Path(input.gff).resolve(),
    container:
        "docker://pasapipeline/pasapipeline:2.5.3"
    shell:
        """
        cd {params.pasa_dir}
        /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta {params.abs_fa} --pasa_transcripts_gff3 {params.abs_gff}
        """


rule cat_all_gff:
    input:
        gff=rules.pasa.output.assemblies_gff,
        blat_gff=rules.pasa.output.blat_gff,
        gmap_gff=rules.pasa.output.gmap_gff,
        transdecoder_gff=rules.pasa_dbi.output.pasa_dbi,
    output:
        all_gff="results/pasa/{species}/{species}.gff3",
    params:
        pasa_dir="results/pasa/{species}",
    shell:
        """
        cat {input.gff} {input.blat_gff} {input.gmap_gff} {input.transdecoder_gff} > {output.all_gff}
        """

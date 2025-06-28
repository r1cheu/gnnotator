
rule run_pasa:
    input:
        pasa_config="results/trinity/{species}/PASA/pasa.alignAssembly.txt",
        transcripts="results/trinity/{species}/PASA/{species}_transcripts.fasta",
        masked_genome="results/mask_repeat/{species}.masked.fa",
        trans_gtf="results/map_rna/{species}/Merge.gtf",
    output:
        pasa_db="results/trinity/{species}/PASA/{species}_DB_pasa.sqlite",
    params:
        pasa_dir="results/trinity/{species}/PASA",
        db_name="{species}_DB_pasa",
    conda:
        "../envs/pasa.yml"
    threads: 32
    shell:
        """
        # PASACONF is a user-defined variable in the original script, pointing to a global config file.
        # PASA should be able to find its config file if it's in the environment.
        # The docker command is replaced by a direct call to the PASA pipeline.
        # The original script uses a hardcoded path to the PASA config, which is not ideal.
        # We will assume PASA is configured correctly in the conda environment.
        Launch_PASA_pipeline.pl -c {input.pasa_config} -C -R -g {input.masked_genome} -t {input.transcripts} --trans_gtf {input.trans_gtf} --ALIGNERS blat,gmap --CPU {threads}
        """

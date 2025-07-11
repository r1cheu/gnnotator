rule create_evm_weight:
    output:
        weight="results/evm/{species}/weight.txt",
    params:
        conf=config["evm_weights"],
    conda:
        "../envs/base.yml"
    log:
        "logs/evm/{species}/create_evm_weight.log",
    script:
        "../scripts/create_evm_weight.py"


rule evidencemodeler:
    input:
        fa=rules.mask_repeat.output.masked,
        transcripts_gff=rules.cat_all_gff.output.all_gff,
        weight=rules.create_evm_weight.output.weight,
        prediction_gff=rules.abinitio_prediction_all.output.all,
        protein_align=rules.miniprot_forEVM.output.evm_gff,
    output:
        cds="results/evm/{species}/{species}.EVM.cds",
    params:
        segment_size=500000,
        overlap_size=50000,
        dir=lambda w: Path(f"results/evm/{w.species}").resolve(),
        fa=lambda w, input: Path(input.fa).resolve(),
        transcripts_gff=lambda w, input: Path(input.transcripts_gff).resolve(),
        weight=lambda w, input: Path(input.weight).resolve(),
        prediction_gff=lambda w, input: Path(input.prediction_gff).resolve(),
        protein_align=lambda w, input: Path(input.protein_align).resolve(),
        log=lambda w: Path(f"logs/evm/{w.species}/evidencemodeler.log").resolve(),
    log:
        "logs/evm/{species}/evidencemodeler.log",
    container:
        container_image("brianjohnhaas/evidencemodeler:2.1.0")
    shell:
        """
        cd {params.dir}
        EVidenceModeler --sample_id {wildcards.species} --genome {params.fa} --weight {params.weight} --gene_predictions {params.prediction_gff} --protein_alignments {params.protein_align} --transcript_alignments {params.transcripts_gff} --segmentSize {params.segment_size} --overlapSize {params.overlap_size} &> {params.log}
        """


rule transeq:
    input:
        cds=rules.evidencemodeler.output.cds,
    output:
        fa="results/evm/{species}/{species}.EVM.aa.fasta",
    log:
        "logs/evm/{species}/transeq.log",
    conda:
        "../envs/emboss.yml"
    shell:
        """
        transeq -sequence {input.cds} -outseq {output.fa} -table 1 &>{log}
        """

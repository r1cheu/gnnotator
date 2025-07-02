rule prepare_evm_weights:
    input:
        weight_template=config["gene_prediction"]["evm_weight_template"],
    output:
        weight_file="results/gene_prediction/{species}/weight.txt",
    shell:
        """
        cp {input.weight_template} {output.weight_file}
        sed -i "s/Nipv2_DB_pasa/{wildcards.species}_DB_pasa/g" {output.weight_file}
        """


rule run_evm:
    input:
        genome="results/mask_repeat/{species}.masked.fa",
        weights="results/gene_prediction/{species}/weight.txt",
        gene_predictions="results/gene_prediction/{species}/abinitio_prediction_all.gff3",
        protein_alignments="results/gene_prediction/{species}/miniprot.gff3",  # Assuming miniprot output is directly usable
        transcript_alignments="results/gene_prediction/{species}/PASA/transcripts.gff3",
    output:
        evm_gff="results/gene_prediction/{species}/{species}.EVM.gff3",
        evm_cds="results/gene_prediction/{species}/{species}.EVM.cds",
        evm_proteins="results/gene_prediction/{species}/{species}.EVM.proteins",
    params:
        sample_id="{species}",
        segment_size=500000,
        overlap_size=50000,
    conda:
        "../envs/gene_prediction.yml"
    threads: 1
    shell:
        """
        EVidenceModeler --sample_id {params.sample_id} \
            --genome {input.genome} \
            --weights {input.weights} \
            --gene_predictions {input.gene_predictions} \
            --protein_alignments {input.protein_alignments} \
            --transcript_alignments {input.transcript_alignments} \
            --segmentSize {params.segment_size} \
            --overlapSize {params.overlap_size} \
            --output_file {output.evm_gff}
        """


rule translate_cds:
    input:
        cds="results/gene_prediction/{species}/{species}.EVM.cds",
    output:
        aa_fasta="results/gene_prediction/{species}/{species}.EVM.aa.fasta",
    conda:
        "../envs/gene_prediction.yml"
    shell:
        """
        transeq -sequence {input.cds} -outseq {output.aa_fasta} -table 1
        """

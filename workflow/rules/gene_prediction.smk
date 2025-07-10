from pathlib import Path


rule fgenesh:
    input:
        fa=rules.mask_repeat.output.masked,
    output:
        fg="results/gene_prediction/{species}/{species}.fgenesh",
    container:
        "docker://leizzzz/ppd:3.1.1"
    log:
        "logs/gene_prediction/{species}/fgenesh.log",
    shell:
        """
        perl /app/run_FgeneSH.pl {input.fa} {output.fg} seq &>{log}
        """


rule fgenesh_name:
    input:
        fg="results/gene_prediction/{species}/{species}.fgenesh",
    output:
        name="results/gene_prediction/{species}/{species}.fgenesh.name",
    container:
        "docker://leizzzz/ppd:3.1.1"
    log:
        "logs/gene_prediction/{species}/fgenesh_name.log",
    shell:
        """
        perl /app/fgenesh_seq_name.pl {input.fg} &>{log}
        """


rule fgenesh_gff3:
    input:
        name="results/gene_prediction/{species}/{species}.fgenesh.name",
    output:
        gff3="results/gene_prediction/{species}/{species}.fgenesh.forEVM.gff3",
    container:
        "docker://brianjohnhaas/evidencemodeler"
    log:
        "logs/gene_prediction/{species}/fgenesh_gff3.log",
    shell:
        """
        perl /usr/local/bin/EvmUtils/misc/fgenesh_to_GFF3.pl {input.name} > {output.gff3} 2>{log}
        """


rule snap:
    input:
        fa=rules.mask_repeat.output.masked,
    output:
        snap="results/gene_prediction/{species}/{species}.snap.gff3",
    container:
        "docker://leizzzz/snap:4ad1e95"
    log:
        "logs/gene_prediction/{species}/snap.log",
    shell:
        """
        snap -gff /app/HMM/O.sativa.hmm {input.fa} > {output.snap} 2>{log}
        """


rule snap_gff3:
    input:
        snap="results/gene_prediction/{species}/{species}.snap.gff3",
    output:
        gff3="results/gene_prediction/{species}/{species}.snap.forEVM.gff3",
    container:
        "docker://brianjohnhaas/evidencemodeler"
    log:
        "logs/gene_prediction/{species}/snap_gff3.log",
    shell:
        """
        perl /usr/local/bin/EvmUtils/misc/SNAP_ExonEtermEinitEsngl_gff_to_gff3.pl {input.snap} > {output.gff3} 2>{log}
        """


rule eviann_prepare:
    input:
        bam="results/map_rna/{species}/Merge.bam",
    output:
        txt="results/gene_prediction/{species}/bam.txt",
    conda:
        "../envs/base.yml"
    threads: 1
    params:
        abs_bam=lambda w, input: Path(input.bam).resolve(),
    log:
        "logs/gene_prediction/{species}/eviann_prepare.log",
    shell:
        """
        echo {params.abs_bam} bam > {output.txt}
        """


rule eviann:
    input:
        fa=rules.mask_repeat.output.masked,
        txt=rules.eviann_prepare.output.txt,
        library="resources/oryza.fasta",
    output:
        gff3="results/gene_prediction/{species}/{species}.masked.fa.pseudo_label.gff",
    container:
        "docker://leizzzz/eviann:2.0.2"
    threads: 20
    log:
        "logs/gene_prediction/{species}/eviann.log",
    params:
        fa=lambda w, input: Path(input.fa).resolve(),
        txt=lambda w, input: Path(input.txt).resolve(),
        library=lambda w, input: Path(input.library).resolve(),
        log=lambda w: Path(f"logs/gene_prediction/{w.species}/eviann.log").resolve(),
    shell:
        """
        cd results/gene_prediction/{wildcards.species}
        eviann.sh -t {threads} -g {params.fa} -r {params.txt} -p {params.library} --liftover &>{params.log}
        """


rule miniprot_index:
    input:
        fa=rules.mask_repeat.output.masked,
    output:
        mpi="results/gene_prediction/{species}/{species}.mpi",
    conda:
        "../envs/miniprot.yml"
    threads: 16
    log:
        "logs/gene_prediction/{species}/miniprot_index.log",
    shell:
        """
        miniprot -t{threads} -d {output.mpi} {input.fa} &> {log}
        """


rule miniprot_gff:
    input:
        mpi=rules.miniprot_index.output.mpi,
        library="resources/IRGSP-1.0.2019-08-29_msu7.0_all_line.fa",
    output:
        gff="results/gene_prediction/{species}/{species}.mpt.gff3",
    conda:
        "../envs/miniprot.yml"
    log:
        "logs/gene_prediction/{species}/miniprot_gff.log",
    threads: 16
    shell:
        """
        miniprot -It{threads} --gff {input.mpi} {input.library} > {output.gff} 2>{log}
        """


rule miniprot_forEVM:
    input:
        gff="results/gene_prediction/{species}/{species}.mpt.gff3",
    output:
        evm_gff="results/gene_prediction/{species}/{species}.mpt.forEVM.gff3",
    container:
        "docker://brianjohnhaas/evidencemodeler"
    log:
        "logs/gene_prediction/{species}/miniprot_forEVM.log",
    shell:
        """
        python3 /usr/local/bin/EvmUtils/misc/miniprot_GFF_2_EVM_GFF3.py {input.gff} > {output.evm_gff} 2>{log}
        """


rule abinitio_prediction_all:
    input:
        fgenesh=rules.fgenesh_gff3.output.gff3,
        snap=rules.snap_gff3.output.gff3,
        eviann=rules.eviann.output.gff3,
    output:
        all="results/gene_prediction/{species}/abinitio_prediction_all.gff3",
    log:
        "logs/gene_prediction/{species}/abinitio_prediction_all.log",
    conda:
        "../envs/base.yml"
    shell:
        """
        cat {input.fgenesh} {input.snap} {input.eviann} | grep -v "^#" > {output.all}
        """

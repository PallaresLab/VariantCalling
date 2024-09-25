rule fastqc:
    input:
        get_fastq
    output:
        zip=config["working_dir"] + "/fastqc/{sample}_RNA_fastqc.zip"
    log:
        log_dir + "/fastqc/{sample}_RNA.log"
    params:
        outdir=config["working_dir"] + "/fastqc/{sample}",
        #old_name=config["working_dir"] + "/fastqc/{sample}_fastqc.zip"
    threads: 32
    conda:
        "../envs/fastqc.yml"
    shell:
        "mkdir -p {params.outdir:q} && "
        "fastqc --quiet "
        "--threads {threads} "
        "--outdir {params.outdir} "
        "--noextract "
        "{config[params][fastqc][extra]} "
        "{input:q} "
        ">{log:q} 2>&1 && "
        "mv '{params.outdir}/'*_fastqc.zip {output:q} && "
        "rm -rf {params.outdir:q} "
    
rule count:
    input:
        gtf=config["ref"]["annotation"],
        bam=config["working_dir"] + "/nudup/{sample}.sorted.markdup.bam" if config["umi"]["mark_umi_duplicates"] else config["working_dir"] + "/star/{sample}/{sample}_RNA_Aligned.sortedByCoord.out.bam"
    output:
        counts=config["working_dir"] + "/featurecounts/{sample}_RNA_counts.tsv",
        summary=config["working_dir"] + "/featurecounts/{sample}_RNA_counts.tsv.summary"
    params:
        ends="" if if_SE  else "-p"
    log:
        log_dir + "/featurecounts/{sample}_RNA.log"
    threads:
        32
    conda:
        "../envs/subread.yml"
    shell:
        "ln -s {input.bam:q} '{wildcards.sample}.bam' && "
        "featureCounts "
        "{params.ends} "
        "-t {config[params][featurecounts][feature_type]} "
        "-g {config[params][featurecounts][attribute_type]} "
        "-a {input.gtf:q} "
        "{config[params][featurecounts][extra]} "
        "-T {threads} "
        "-o {output.counts:q} "
        "'{wildcards.sample}.bam' "
        ">{log:q} 2>&1 && "
        "rm '{wildcards.sample}.bam'"

rule combined_counts:
    input:
        expand(config["working_dir"] + "/featurecounts/{sample}_RNA_counts.tsv",
               sample = SAMPLES_ALL if if_SE else SAMPLES_PAIRED)
    output:
        config["results_dir"] + "/combined_gene_counts.tsv"
    run:
        shell("cut -f 1,7 {input[0]:q} > {output:q}")
        for f in input:
            shell("cp {output:q} tmp")
            shell("cut -f 7 {f:q} | paste tmp - > {output:q}")
        shell("rm tmp")

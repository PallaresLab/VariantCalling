rule samtools:
    input:
        config['reference_genome']

    output:
        fai = working_dir+"/genome_prepare/"+ref_basename+".fai",
        new_genome = working_dir+"/genome_prepare/"+ref_basename
        
    log:
        log_dir+"/genome_prepare/samtools_faidx.log"
        
    params:
        outdir = working_dir+"/genome_prepare",
        new_genome = working_dir+"/genome_prepare/"+ref_basename 
    resources:
        mem_mb = lambda wildcards, attempt: 30000 * (2 ** (attempt - 1))
        
    conda:
        "../envs/samtools.yml"

    shell: 
        "mkdir -p {params.outdir} && "
        "cp {input} {params.outdir} && "
        "samtools faidx {params.new_genome} "
        ">{log} 2>&1"
        
rule CreateSequenceDictionary:
    input:
        R = working_dir+"/genome_prepare/"+ref_basename,
        dbSNP=config['dbSNP']

    output:
        working_dir+"/genome_prepare/"+ref_name +".dict",
           
    log:
        log_dir+"/genome_prepare/gatk_CreateSequenceDictionary.log"
    resources:
        mem_mb = lambda wildcards, attempt: 30000 * (2 ** (attempt - 1))
        
    conda:
        "../envs/gatk.yml"

    shell: 
        "gatk CreateSequenceDictionary -R {input.R} -O {output} "
        ">{log} 2>&1 && "
        "gatk IndexFeatureFile -I {input.dbSNP} "
        ">>{log} 2>&1"
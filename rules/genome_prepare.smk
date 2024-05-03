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
        
    conda:
        "../envs/samtools.yml"  

    shell: 
        "mkdir -p {params.outdir} && "
        "cp {input} {params.outdir} && "
        "samtools faidx {params.new_genome} "
        ">{log} 2>&1"
        
rule CreateSequenceDictionary:
    input:
        working_dir+"/genome_prepare/"+ref_basename

    output:
        working_dir+"/genome_prepare/"+ref_name +".dict",
           
    log:
        log_dir+"/genome_prepare/gatk_CreateSequenceDictionary.log"
                
    conda:
        "../envs/gatk.yml"

    shell: 
        "gatk CreateSequenceDictionary -R {input} -O {output} "
        ">{log} 2>&1"
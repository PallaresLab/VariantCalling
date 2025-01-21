rule GenomicsDBImport:
    input:
        expand(working_dir+"/HaplotypeCaller/{sample}/{interval}.g.vcf.gz", sample=samples.keys(), interval="{interval}"),
    
    output:
        directory("/dev/shm/pallares_lab/database_{interval}"),
    
    log:
        log_dir + "/JointCallSNPs/GenomicsDBImport_{interval}.log",   
    
    params:
        outdir = "/dev/shm/pallares_lab",
        interval = "{interval}",
        input_params = lambda wildcards, input: " ".join(["-V "+ f for f in input]),
    
    threads: 4
    
    conda:
        "../envs/gatk.yml"
    
    shell:
        "mkdir -p {params.outdir} && "
        "gatk --java-options '-Xmx100G' "
        "GenomicsDBImport "
        "--genomicsdb-workspace-path {output} "
        "-L {params.interval} "
        "{params.input_params} "
        "--batch-size 50 "
        "--reader-threads {threads} "
        ">{log} 2>&1 "

rule GenotypeGVCFs:
    input:
        R = working_dir + "/genome_prepare/"+ref_basename,
        db = "/dev/shm/pallares_lab/database_{interval}",
    
    output:
        working_dir + "/JointCallSNPs/{interval}.vcf.gz",
    
    log:
        log_dir + "/JointCallSNPs/GenotypeGVCFs_{interval}.log",     
        
    threads :1
    
    conda:
        "../envs/gatk.yml"
    
    shell:     
        "gatk --java-options '-Xmx100G' "
        "GenotypeGVCFs "
        "-R {input.R} "
        "-V gendb://{input.db} "
        "-O {output} && "
        "rm -rf {input.db} "
        ">{log} 3>&1 "


rule GatherVcfs:
    input:
       expand(working_dir + "/JointCallSNPs/{interval}.vcf.gz", interval=INTERVALS)
    output:
        working_dir + "/JointCallSNPs/all_samples.vcf.gz"
    log:
        log_dir + "/JointCallSNPs/GatherVcfs.log",     
    
    params:
        input_params = lambda wildcards,input:" ".join(["I= "+ f for f in input])
        
    conda:
        "../envs/gatk.yml"

    shell: 
        "rm -rf /dev/shm/pallares_lab/ && "    
        "gatk GatherVcfs "
        "{params.input_params} O={output} && "
        "rm {input} && "
        "gatk IndexFeatureFile -I {output} "
        ">{log} 3>&1 " 
     
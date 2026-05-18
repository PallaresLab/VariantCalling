rule GenomicsDBImport:
    input:
        expand(working_dir+"/HaplotypeCaller/{sample}/{{idx}}.g.vcf.gz", sample=samples.keys()),

    
    output:
        directory(working_dir +"/JointCallSNPs/database_{idx}"),
    
    log:
        log_dir + "/JointCallSNPs/GenomicsDBImport_{idx}.log",   
    
    params:
        outdir = working_dir +"/JointCallSNPs",
        input_params = lambda wildcards, input: " ".join(["-V "+ i for i in input]),
        interval = lambda wildcards: " ".join([f"-L {chrom}:{start+1}-{end}" for chrom, start, end in INTERVALS[int(wildcards.idx)]]),
        
    threads: 4

    resources:
        mem_mb = lambda wildcards, attempt: 20000 * (2 ** (attempt - 1))

    conda:
        "../envs/gatk.yml"

    shell:
        "mkdir -p {params.outdir} && "
        "gatk --java-options '-Xmx8G' "
        "GenomicsDBImport "
        "--genomicsdb-workspace-path {output} "
        "{params.interval} "
        "{params.input_params} "
        "--batch-size 50 "
        "--reader-threads {threads} "
        ">{log} 2>&1 "

rule GenotypeGVCFs:
    input:
        R = working_dir + "/genome_prepare/"+ref_basename,
        db = working_dir +"/JointCallSNPs/database_{idx}",
    
    output:
        working_dir + "/JointCallSNPs/{idx}.vcf.gz",
    
    log:
        log_dir + "/JointCallSNPs/GenotypeGVCFs_{idx}.log",     
        
    threads: 1

    resources:
        mem_mb = lambda wildcards, attempt: 20000 * (2 ** (attempt - 1))

    conda:
        "../envs/gatk.yml"

    shell:
        "gatk --java-options '-Xmx8G' "
        "GenotypeGVCFs "
        "-R {input.R} "
        "-V gendb://{input.db} "
        "-O {output} "
        ">{log} 2>&1 && "
        "rm -rf {input.db} "


rule GatherVcfs:
    input:
       expand(working_dir + "/JointCallSNPs/{idx}.vcf.gz", idx=IDX)
    output:
        working_dir + "/JointCallSNPs/all_samples.vcf.gz"
    log:
        log_dir + "/JointCallSNPs/GatherVcfs.log",     
    
    params:
        input_params = lambda wildcards,input:" ".join(["I= "+ f for f in input])
        
    resources:
        mem_mb = lambda wildcards, attempt: 50000 * (2 ** (attempt - 1))

    conda:
        "../envs/gatk.yml"

    shell:
        "gatk GatherVcfs "
        "{params.input_params} O={output} "
        ">{log} 2>&1 && "
        "rm {input} && "
        "gatk IndexFeatureFile -I {output} "
        ">>{log} 2>&1 "
     
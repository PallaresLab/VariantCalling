rule HaplotypeCaller:
    input:
        R = working_dir+"/genome_prepare/"+ref_basename,
        table = working_dir+"/BQSR/{sample}.table",
        recal_bam = working_dir+"/BQSR/{sample}_recal.bam",
    
    output:
        working_dir+"/HaplotypeCaller/{sample}.g.vcf.gz",
        
    
    log:
        log_dir+"/HaplotypeCaller/{sample}.log",

        
    params:
        outdir=working_dir+"/HaplotypeCaller"
    
        
    shell:
        "mkdir -p {params.outdir} && "
        "gatk --java-options '-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4' "
        "HaplotypeCaller "
        "-O {output} "
        "-R {input.R} "
        "-I {input.recal_bam} "
        "-ERC GVCF "      
        ">{log} 2>&1 "

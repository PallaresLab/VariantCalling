rule HaplotypeCaller:
    input:
        R = working_dir+"/genome_prepare/"+ref_basename,
        table = working_dir+"/BQSR/{sample}.table",
        recal_bam = working_dir+"/BQSR/{sample}_recal.bam",
    
    output:
        working_dir+"/HaplotypeCaller/{sample}/{interval}.g.vcf.gz",
        
    
    log:
        log_dir+"/HaplotypeCaller/{sample}/{interval}.log",

        
    params:
        outdir=working_dir+"/HaplotypeCaller/{sample}",
        interval ="{interval}",
    
    threads :1
    
    conda:
        "../envs/gatk.yml"
  
    shell:
        "mkdir -p {params.outdir} && "
        "gatk --java-options '-Xmx20G' "
        "HaplotypeCaller "
        "-O {output} "
        "-R {input.R} "
        "-I {input.recal_bam} "
        "-L {params.interval} "
        "-ERC GVCF "      
        ">{log} 2>&1 "
        
        
        
   

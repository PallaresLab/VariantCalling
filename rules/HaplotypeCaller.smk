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
    
    threads: multiprocessing.cpu_count()
        
    conda:
        "../envs/gatk.yml"
        
    shell:
        "mkdir -p {params.outdir} && "
        "gatk --java-options '-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4' "
        "HaplotypeCallerSpark "
        "-O {output} "
        "-R {input.R} "
        "-I {input.recal_bam} "
        "-ERC GVCF "
        "--spark-master local[{threads}] "        
        ">{log} 2>&1 "
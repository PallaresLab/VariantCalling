rule VariantFiltration:
    input:
        V = working_dir + "/VQSR/SNP_recal.vcf.gz",
    
    output:
        working_dir + "/HardFilter/all_samples_filter.vcf.gz",

    log:
        log_dir+"/HardFilter/VariantFiltration.log",

        
    params:
        outdir=working_dir+"/HardFilter"
        
    conda:
        "../envs/gatk.yml"
        
    shell:
        "mkdir -p {params.outdir} && "
        "gatk --java-options '-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4' "
        "VariantFiltration "
        "-O {output} "
        "-V {input.V} "
        "--filter-name 'DPfilter' "
        "--filter-expression 'DP<10' "      
        ">{log} 2>&1 "
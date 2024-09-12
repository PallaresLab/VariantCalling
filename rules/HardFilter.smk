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
        "--filter-name 'LowQD' --filter-expression 'QD < 5.0' "
        "--filter-name 'LowMQ' --filter-expression 'MQ < 50.0' "
        "--filter-name 'HighStrandBiasFS' --filter-expression 'FS > 10.0' "
        "--filter-name 'HighStrandBiasSOR' --filter-expression 'SOR > 2.0' "
        "--filter-name 'ExtremeDepth' --filter-expression 'DP < 10 || DP > 200' "
        "--filter-name 'LowPE' --filter-expression 'PE < 500000' "
        ">{log} 2>&1 "

        #Initial Hard Filtering on the VCF:
        #Quality by Depth (QD): Filter out SNPs with low QD (e.g., QD < 3 or 5) to remove variants that are likely low confidence
        #Mapping Quality (MQ): Apply a filter to remove SNPs with low MQ (e.g., MQ < 50) to exclude poorly mapped reads
        #Strand Bias (FS and SOR): Filter out SNPs with significant strand bias (e.g., FS > 10 or SOR > 2) to reduce false positives caused by strand-specific errors.
        #Depth of Coverage (DP): Remove SNPs with extreme coverage (too high or too low), as these may represent mapping artifacts or low-confidence calls.
        #the last parameter lowPE discards anything that is lower than 500k reads per sample.


rule SelectVariants:
    input:
        filtered_vcf = working_dir + "/HardFilter/all_samples_filter.vcf.gz",

    output:
        selectVariants_vcf = working_dir + "/HardFilter/all_samples_selectVariants.vcf.gz",

    log:
        log_dir + "/HardFilter/SelectVariants.log",

    params:
        outdir = working_dir + "/HardFilter"

    conda:
        "../envs/gatk.yml"

    shell:
        "mkdir -p {params.outdir} && "
        "gatk --java-options '-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4' "
        "SelectVariants "
        "-V {input.filtered_vcf} "
        "--set-filtered-genotype-to-no-call "
        "--exclude-filtered "
        "--genotype-filter-expression 'GQ < 9' "
        "--genotype-filter-name 'LowGQ' "
        "-O {output.selectVariants_vcf} "
        ">{log} 2>&1 "

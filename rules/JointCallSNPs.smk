rule GenomicsDBImport:
    input:
        f = expand(working_dir+"/HaplotypeCaller/{sample}.g.vcf.gz", sample=samples.keys()),
        L = config['bed_file'] 
    
    output:
        directory(working_dir + "/JointCallSNPs/database"),
    
    log:
        log_dir + "/JointCallSNPs/GenomicsDBImport.log",   
     
    
    params:
        outdir = working_dir+"/JointCallSNPs"
        
    threads:
        32
        
    conda:
        "../envs/gatk.yml"

    shell:
        "mkdir -p {params.outdir} && "
        "echo {input.f} | sed 's/ /\\n /g' | sed 's/^/-V /' > input_files.txt && "
        
        "gatk --java-options '-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4' "
        "GenomicsDBImport "
        "--genomicsdb-workspace-path {output} -L {input.L} "
        "$(< input_files.txt) " #info about -V *.g.vcf.gz
        "--reader-threads {threads} "
        "--max-num-intervals-to-import-in-parallel 4 "
        ">{log} 2>&1 && "
        
        "rm input_files.txt  "
        
rule GenotypeGVCFs:
    input:
        R = working_dir + "/genome_prepare/"+ref_basename,
        db = working_dir + "/JointCallSNPs/database",
    
    output:
        working_dir + "/JointCallSNPs/all_samples.vcf.gz",
    
    log:
        log_dir + "/JointCallSNPs/GenotypeGVCFs.log",     
    
    params:
        outdir = working_dir+"/JointCallSNPs"
        
    conda:
        "../envs/gatk.yml"

    shell:     
        "gatk --java-options '-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4' "
        "GenotypeGVCFs "
        "-R {input.R} "
        "-V gendb://{input.db} "
        "-O {output} "
        ">{log} 2>&1"
        
        




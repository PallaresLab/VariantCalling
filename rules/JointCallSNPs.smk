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
        256
        
    conda:
        "../envs/gatk.yml"

    shell:
        "mkdir -p {params.outdir} && "
        "echo {input.f} | sed 's/ /\\n /g' | sed 's/^/-V /' > input_files.txt && "
        
        "gatk --java-options '-Xmx100G -XX:+UseParallelGC -XX:ParallelGCThreads=32' "
        "GenomicsDBImport "
        "--genomicsdb-workspace-path {output} "
        "-L out.bed "
        "$(< input_files.txt) " #info about -V *.g.vcf.gz
        #"--reader-threads {threads} "
        "--max-num-intervals-to-import-in-parallel 32 "
        "--batch-size 50 "
        "--reader-threads 5 "
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
    threads:
        256

    shell:     
        "gatk --java-options '-Xmx100G -XX:+UseParallelGC -XX:ParallelGCThreads=32' "
        "GenotypeGVCFs "
        "-R {input.R} "
        "-V gendb://{input.db} "
        "-O {output} "
        ">{log} 2>&1"
        
        




shell.executable("bash")

from snakemake.utils import min_version
#min_version("8.10.6")

configfile: "config.yaml"
log_dir = config['output_dir'] + "/logs"
working_dir = config['output_dir'] + "/working"
ref_genome_name=os.path.basename(config['reference_genome'])
sample_files = snakemake.utils.listfiles(config["bam_dir"]+"/{sample}.bam")
samples = dict((y[0], x) for x, y in sample_files)




rule all:
    input:
        directory(working_dir + "/JointCallSNPs/database"),
        #working_dir + "/JointCallSNPs/all_samples.vcf.gz",
     
'''
include: "rules/genome_prepare.smk"         

        
        
rule bqsr:
    input:
        bam = config['bam_dir']+"/{sample}.bam",
        R = working_dir + "/genome/"+ref_genome_name,
        vcf = config['dbSNP']
    
    output:
        table=working_dir+"/bqsr/{sample}.table",
        recal_bam=working_dir+"/bqsr/{sample}_recal.bam"
        
    log:
        log1 = log_dir+"/bqsr/BaseRecalibrator/{sample}.log",
        log2 = log_dir+"/bqsr/ApplyBQSR/{sample}.log"
        
        
    params:
        outdir = working_dir+"/bqsr"
        
    threads: 32
        
    #conda:
        #"envs/gatk.yml"

    shell:
        "mkdir -p {params.outdir} &&"
        "gatk BaseRecalibratorSpark "
        "-O {output.table} "
        "-R {input.R} "
        "-I {input.bam} "
        "--known-sites {input.vcf} "
        "--spark-master local[{threads}] "        
        ">{log.log1} 2>&1 &&"
        
        "gatk ApplyBQSRSpark "
        "-O {output.recal_bam} "
        "-R {input.R} "
        "-I {input.bam} "
        "--bqsr-recal-file {output.table} "
        "--spark-master local[{threads}] "
        ">{log.log2} 2>&1"
  

rule HaplotypeCaller:
    input:
        R=working_dir+"/genome/"+ref_genome_name,
        table=working_dir+"/bqsr/{sample}.table",
        recal_bam=working_dir+"/bqsr/{sample}_recal.bam",
    
    output:
        gvcf = working_dir+"/HaplotypeCaller/{sample}.g.vcf.gz",
        
    
    log:
        log_dir+"/HaplotypeCaller/{sample}.log",

        
    params:
        outdir=working_dir+"/HaplotypeCaller"
    
    threads: 32
        
    #conda:
        #"envs/gatk.yml"

    shell:
        "mkdir -p {params.outdir} && "
        "gatk HaplotypeCallerSpark "
        "-O {output.gvcf} "
        "-R {input.R} "
        "-I {input.recal_bam} "
        "-ERC GVCF "
        "--spark-master local[{threads}] "        
        ">{log} 2>&1 "


'''
rule JointCallSNPs:
    input:
        f = expand(working_dir+"/HaplotypeCaller/{sample}.g.vcf.gz", sample=samples.keys()),
        R = working_dir + "/genome/"+ref_genome_name,
        L = config['bed_file'] 
    
    output:
        db = directory(working_dir + "/JointCallSNPs/database"),
        #vcf = working_dir + "/JointCallSNPs/all_samples.vcf.gz",
    
    log:
        log1 = log_dir + "/JointCallSNPs/GenomicsDBImport.log",   
        #log2 = log_dir + "/JointCallSNPs/GenotypeGVCFs.log",     
    
    params:
        outdir = working_dir+"/JointCallSNPs"
        
    #conda:
        #"envs/gatk.yml"

    shell:
        "mkdir -p {params.outdir} && "
        "echo {input.f} | sed 's/ /\\n /g' | sed 's/^/-V /' > input_files.txt && "
        
        "gatk GenomicsDBImport "
        "--genomicsdb-workspace-path {output.db} -L {input.L} "
        "$(< input_files.txt) "#info about -V *.g.vcf.gz
        ">{log.log1} 2>&1 && "
        
        "rm input_files.txt  "
'''       
        "gatk GenotypeGVCFs "
        "-R {input.R} "
        "-V gendb://{output.db} "
        "-O {output.vcf} "
        ">{log.log2} 2>&1"
        
        
        
rule VQSR:
    input:
        V = working_dir + "/JointCallSNPs/all_samples.vcf.gz",
        R = working_dir + "/genome/"+ref_genome_name,
        L = config['bed_file'] 
    
    output:
        output1 = working_dir + "/VQSR/SNP.recal",
        output2 = working_dir + "/VQSR/SNP_recal.vcf.gz",

    
    log:
        log1 = log_dir + "/VQSR/VariantRecalibrator.log",
        log2 = log_dir + "/VQSR/ApplyVQSR.log",
   
    
    params:
        outdir = working_dir+"/VQSR"
        
    #conda:
        #"envs/gatk.yml"

    shell:
        "mkdir -p {params.outdir} && "
        "gatk VariantRecalibrator "
        "-V {input.V} "
        "-R {input.R} "
        "-L {input.L} "
        "-resource:dgrp2,known=false,training=true,truth=true,prior=15.0 dgrp2.vcf "
        "-resource:gdl,known=false,training=true,truth=true,prior=15.0 gdl.vcf "
        "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR "
        "-mode SNP "
        "--tranches-file output.tranches "
        " --rscript-file recalibration_plots.R "
        "-O {output.output1} "
        ">{log.log1} 2>&1 && "
        
        "gatk ApplyVQSR "
        "-V {input.V} "
        "-R {input.R} "
        "--tranches-file output.tranches "
        "-mode SNP "
        "--recal-file SNP.recal "
        "--truth-sensitivity-filter-level 99.5"
        "-O {output.output2} "
        ">{log.log2} 2>&1"
'''        
       
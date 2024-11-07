shell.executable("bash")

from snakemake.utils import min_version
import multiprocessing
min_version("8.10.7")

configfile: "config.yaml"
log_dir = config['output_dir'] + "/logs"
working_dir = config['output_dir'] + "/working"
ref_basename=os.path.basename(config['reference_genome'])
ref_name=os.path.splitext(os.path.basename(config['reference_genome']))[0]    
sample_files = snakemake.utils.listfiles(config["bam_dir"]+"/{sample}.bam")
samples = dict((y[0], x) for x, y in sample_files)

singularity: "docker://broadinstitute/gatk"


rule all:
    input:
        working_dir + "/HardFilter/all_samples_filter.vcf.gz",
        
        
#include: "rules/genome_prepare.smk"        
#include: "rules/BQSR.smk"
#include: "rules/HaplotypeCaller.smk"        
include: "rules/JointCallSNPs.smk"
include: "rules/VQSR.smk"
include: "rules/HardFilter.smk"
     

       
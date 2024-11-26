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

def split_range(n, max_range):
    ranges = []
    for start in range(0, n, max_range):
        end = min(start + max_range, n)
        ranges.append((start, end))
    return ranges
def split_names_by_length(names_dict, max_length):
    result = {}
    for name, length in names_dict.items():
        result[name] = split_range(length,max_length)
    return result
    
genomic_dict={}
    
with open(config['bed_file']) as f:
    for line in f:
        spt=line.strip().split("\t")
        if len(spt)==3:
            genomic_dict[spt[0]]=int(spt[2])
gd=split_names_by_length(genomic_dict,100000)

INTERVALS = [f"{name}:{start+1}-{end}" for name, length in gd.items() for start, end in length]


rule all:
    input:
        expand(working_dir + "/HaplotypeCaller/{sample}/{interval}.g.vcf.gz", sample=samples.keys(), interval=INTERVALS),
        working_dir + "/VQSR/SNP_recal.vcf.gz",
        
        
include: "rules/genome_prepare.smk"        
include: "rules/BQSR.smk"
include: "rules/HaplotypeCaller.smk"        
include: "rules/JointCallSNPs.smk"
include: "rules/VQSR.smk"

     

     

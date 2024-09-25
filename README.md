# Variant Calling pipeline

## Description
This snakemake pipeline is designed for Variant Calling

### Inputs

* the bam files
* reference genome
* SNP database
* bed file


### Outputs

*   `logs\` - Directory of log files for each job, check here first if you run into errors
*   `working\` - Directory containing intermediate files for each job

### Workflow

1.  **genome_prepare 
2.  **BQSR--BaseRecalibrator, ApplyBQSR
3.  **HaplotypeCaller--HaplotypeCaller
4.  **JointCallSNPs--JointCallSNP, GenotypeGVCFs
5.  **VQSR--VariantRecalibrator, ApplyVQSR
7.  **HardFilter--VariantFiltration


## Setup environment

1.  Install conda

2.  Clone workflow into working directory

    ```bash
    git clone <repo> <dir>
    cd <dir>
    ```
    
3.  Create a new enviroment

    ```bash
    conda env create -n <project> --file environment.yaml
    ```

3.  Activate the environment

    ```bash
    conda activate <project_name>
    ```

4. Install snakemake

    ```bash
    conda install snakemake
    ```

5.  Edit configuration files
    change the path of bam_dir, output_dir, reference_genome, dbSNP and bed_file in "config.yaml"

6.  Create index for your SNP database

    ```bash
    gatk IndexFeatureFile -I dbSNP.vcf.gz
    ```

7.  The first time you are executing this snakemake pipeline it should run locally, once the first run is over (you can use --dry), you can switch to running it on the cluster.

    ```bash
    snakemake --configfile "config.yaml" --use--singularity  --cores N --dryrun
    ```

8.  Execute the workflow

    ```bash
    snakemake --configfile "config.yaml" --use--singularity  --cores N
    ```



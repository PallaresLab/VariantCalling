# Paired-end DNA pipeline

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

4.  Enable the [Bioconda](https://bioconda.github.io/#using-bioconda) channel

    ```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    ```

5. Install snakemake

    ```bash
    conda install snakemake
    ```

6.  Edit configuration files
    change the path of fastq_dir, output_dir, reference_genome in "config.yaml"

7.  Create index for your SNP database

    ```bash
    gatk IndexFeatureFile -I dbSNP.vcf.gz
    ```

8.  The first time you are executing this snakemake pipeline it should run locally, once the first run is over (you can use --dry), you can switch to running it on the cluster.

    ```bash
    snakemake --configfile "config.yaml" --use-conda  --cores N --dryrun
    ```

9.  Execute the workflow

    ```bash
    snakemake --configfile "config.yaml" --use-conda  --cores N
    ```



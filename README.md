

# nf-qc
ABMicroBioinf/nf-qc is a bioinformatics pipeline that can be used to do the genome/metagenome sequence quality control and generate sequence stats information. For now, I only tested the pipeline with the illumina paired-end dataset. I will keep updating the program later.

![flowchart](https://user-images.githubusercontent.com/52679027/121734696-5e4c9a00-cab2-11eb-80cb-fa9730cfaadc.png)

# Third-party software
This pipeline are depending on a number of the third-party software. Please install the following 3rd party dependencies and make sure they are on your system path
* nextflow
* bbmap
* fastqc
* multiqc
* kraken2
* bracken
* seqtk

# Quick Start
To run on the cluster in the interactive mode, you need to load the requried module first. Example module loading from our system:
```
module load nextflow/v20.10.0
module load bbmap/v38.90
module load miniconda3/fastqc
module load local_conda_env/multiqc
module load kraken2/v2.1.1
module load bracken/v2.6
module load seqtk/v1.3-r117-dirty
```

Run the program locally and get help:
```
nextflow run main.nf --help
```

Example command to run the program from github
```
nextflow run ABMicroBioinf/nf-qc --inputdir "fastq" --paired_end "*_{1,2}.fastq"  --minlen "20" -r c40d082
```
In the command, "c40d082" is the github revision number and can change. This command processes the fastq format files contained inside "fastq" directory. The files contained inside the directory are:
* ERR019555_1.fastq
* ERR019555_2.fastq
* ERR019564_1.fastq
* ERR019564_2.fastq

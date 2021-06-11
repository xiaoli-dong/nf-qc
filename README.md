
# nf-qc
ABMicroBioinf/nf-qc is a bioinformatics pipeline that can be used to do the genome/metagenome sequence quality control and generate sequence stats information

![flowchart](https://user-images.githubusercontent.com/52679027/121611912-9a322180-ca16-11eb-9447-663dd6ffd4af.png)

# Third-party software
This pipeline are depending on a number of the third-party software. Please install the following 3rd party dependencies and make sure they are on your system path
* nextflow/v20.10.0
* bbmap/v38.90
* fastqc
* multiqc
* kraken2/v2.1.1
* bracken/v2.6
* seqtk/v1.3-r117-dirty

# Quick Start
To run on the cluster in the interactive mode, you need to load the requried module first:
'''
module load nextflow/v20.10.0
module load bbmap/v38.90
module load miniconda3/fastqc
module load local_conda_env/multiqc
module load kraken2/v2.1.1
module load bracken/v2.6
module load seqtk/v1.3-r117-dirty
'''
Get help:
'''
nextflow run main.nf --help
'''

Example command for processing the fastq format paired-end sequence files(ERR019555_1.fastq  ERR019555_2.fastq  ERR019564_1.fastq  ERR019564_2.fastq) contained inside the fastq directory
'''
nextflow run ABMicroBioinf/nf-qc --inputdir "fastq" --paired_end "*_{1,2}.fastq" -with-dag flowchart.png --minlen "20" -r c40d082

'''
c40d082 is the github revision and can change
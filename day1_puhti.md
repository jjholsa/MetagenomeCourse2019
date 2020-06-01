# Day 1
*Antti Karkman, Katariina Pärnänen and Jenni Hultman*

# Metagenome analysis of infant gut metagenomes
Log into Taito, either with ssh (Mac/Linux) or PuTTy (Windows)

__All the scripts are to be run in `Metagenomics2019` folder!__

## Data download
First set up the course directory, make some folders and then download the data.  
Move to work directory at CSC and make a course directory there. Or a subfolder under a previous course folder.  
```
cd /scratch/project_xxxx/
mkdir Metagenomics2019
cd Metagenomics2019
mkdir raw_data
mkdir scripts
mkdir trimmed_data
```

Download the metagenomic data (takes few minutes)  This is data Igor gave you the path for
```
cd raw_data
cp /wrk/parnanen/shared/COURSE_DATA/* .
```

Make a file containing the sample names. This is the second field separated by `_`.    
Use only the forward reads for this.  
 `ls *_R1.fastq.gz |awk -F "_" '{print $2}' > ../sample_names.txt`  

## QC and trimming
QC for the raw data (takes few min, depending on the allocation).  
Go to the folder that contains the raw data and make a folder called e.g. `FASTQC` for the QC reports.  
Then run the QC for the raw data and combine the results to one report using `multiqc`.  

Can be done on the interactive nodes using `sinteractive`. see https://docs.csc.fi/computing/running/interactive-usage/
Just go to the right folder, activate `QC_env` and run `FastQC`. After it is finished you can just log out from the computing node. You don´t need to free any resources.

```
# start interactive computing node, change xxx to your Puhti project number

sinteractive --account project_xxxxxx --time 04:00:00 --mem 8000 --tmp 100

# activate the QC environment
module load bioconda/3
source activate QC_env

# Run fastqc
fastqc ./*.fastq.gz -o FASTQC/ -t 4

# Then combine the reports with multiqc
multiqc ./ --interactive

# deactivate the virtual env
source deactivate

# log out from the computing node
exit

# and free the resources after the job is done
exit
```

Copy the resulting HTML file to your local machine with `scp` from the command line (Mac/Linux) or *WinSCP* on Windows. Have a look at the QC report with your favourite browser.  

After inspecting the output, it should be clear that we need to do some trimming.  
__What kind of trimming do you think should be done?__

You can check the adapter sequences from Illumina's website. (e.g. search with "*Illumina adapter sequences*"). Or from [here](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-10.pdf).

For trimming we'll make a bash script that runs `cutadapt` for each file using the `sample_names.txt` file.    
Go to your scripts folder and make a bash script for cutadapt with any text editor. Specify the adapter sequences that you want to trim after `-a` and `-A`. What is the difference with `-a` and `-A`? And what is specified with option `-O`? You can find the answers from Cutadapt [manual](http://cutadapt.readthedocs.io).

Save the file as `cutadapt.sh` in the scripts folder.  
__Make sure that the PATHS are correct in your own script__  

```
#!/bin/bash

while read i
do
        cutadapt  -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -O 10 --max-n 0 \
        -o ../trimmed_data/$i"_R1_trimmed.fastq" -p ../trimmed_data/$i"_R2_trimmed.fastq" \
        *$i*_R1.fastq.gz *$i*_R2.fastq.gz > ../trimmed_data/$i"_trim.log"
done < $1
```

Then we need a batch job file to submit the job to the SLURM system. More about CSC batch jobs here: https://research.csc.fi/taito-batch-jobs  
Make another file with text editor that runs the script above.

```
#!/bin/bash -l

#SBATCH --job-name=cutadapt
#SBATCH --account=project_xxxxx
#SBATCH --time=02:00:00
#SBATCH --mem=1000
#SBATCH --partition=small
#SBATCH --output=cutadapt_out_%j.txt
#SBATCH --error=cutadapt_err_%j.txt
#

module load biokit
cd /scratch/project_xxxx/Metagenomics2019/raw_data
bash ../scripts/cutadapt.sh ../sample_names.txt
```

After it is done, we can submit it to the SLURM system. Do it from the course main folder, so go one step back in your folders.  

`sbatch scripts/cut_batch.sh`  

You can check the status of your job with:  

`squeue -l -u $USER`  

After the job has finished, you can see how much resources it actually used and how many billing units were consumed. `JOBID` is the number after the batch job error and output files.  

`seff JOBID`  

Then let's check the results from the trimming.

Go to the folder containing the trimmed reads (`trimmed_data`) and make a new folder (`FASTQC`) for the QC files.  


```
# start interactive computing node, change xxx to your Puhti project number

sinteractive --account project_xxxxxx --time 04:00:00 --mem 8000 --tmp 100

# activate the QC environment
module load bioconda/3
source activate QC_env

# run QC on the trimmed reads
fastqc ./*.fastq -o FASTQC/ -t 4
multiqc ./ --interactive

# deactivate the virtual env
source deactivate

# log out from the computing node
exit

# and free the resources after the job is done
exit
```

Copy it to your local machine as earlier and look how well the trimming went.  

## Assembly
We will assemble all 10 samples together (co-assembly) and use [Megahit assembler](https://github.com/voutcn/megahit) for the job. In addition, we will use MetaQuast to get some statistics about our assembly.  

Megahit is an ultra fast assembly tool for metagenomics data. It is installed to CSC and be loaded with following commands:
```
module purge

module load biokit
```
module purge is needed to remove wrong Python versions you might have loaded earlier today.

Assembling metagenomic data can be very resource demanding and we need to do it as a batch job. As we want to do both individual and co-assemblies the R1 and R2 reads need to be merged into two files with `cat`

```
cd trimmed_data
cat *R1_trimmed.fastq > all_R1.fastq
cat *R2_trimmed.fastq > all_R2.fastq
```
Make a script called co_assembly.sh in a text editor
```
#!/bin/bash
#SBATCH --job-name=megahit
#SBATCH --account=project_xxxxx
#SBATCH --time=24:00:00
#SBATCH --mem=60000
#SBATCH --partition=small
#SBATCH --output=megahit_out_%j.txt
#SBATCH --error=megahit_err_%j.txt
#


module purge

module load biokit

cd /scratch/project_xxxx/Metagenomics2019/

megahit -1 trimmed_data/all_R1.fastq -2 trimmed_data/all_R2.fastq \
         -o co-assembly -t $SLURM_CPUS_PER_TASK --min-contig-len 1000

# MetaQUAST assembly statistics
module purge
module load biokit
cd co-assembly
metaquast.py -t $SLURM_CPUS_PER_TASK --no-plots -o assembly_QC final.contigs.fa
```
Submit the batch job as previously



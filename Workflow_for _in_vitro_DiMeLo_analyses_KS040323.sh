#Workflow for in vitro DiMeLo analyses
#04.03.2023
#Authored by Kousik Sundararajan

### README

### The following code is for analyses of in vitro DiMeLo sequencing data output by ONT in the form of fast5 files. The final output is in the form of bam file for each barcode, with reads filtered by length.
### The code below corresponds to the following workflow

#Input (from nanopore) - fast5 (filename.fast5)
#Fast basecalling (Guppy) - This step is only for demultiplexing the fast5's for faster processing times
#fast5 --> fastq
#Demultiplexing (fastq --> demultiplexed fastq's in different folders, list of readid's per barcode)
#splitting fast5's (using readid's)
#Methylation sensitive basecalling and mapping to template (Megalodon, uses Guppy)
##Corrects sequence based on template
##Maps sequence onto template
##Calls methylation at A's and C's
#fast5 --> fastq and bam files
#bam files have both sequence and methylation information that we will read out using Pysam in a Jupyter notebook

### In the following lines of code, replace <text> with the correct filepaths, filenames, parameters etc.

### All the code below is in bash (or to be pasted in the terminal). The Guppy and Megalodon commands run fastest in the GPU node.
### These commands can be directly pasted into the terminal or run as sbatch scripts. If running as sbatch scripts, do the following:


## 1. First vim the file <filename>.sbatch

vim <sbatch_filename>.sbatch

## 2. Enter the following after pressing I for insert (starting from "#!/bin/bash -l")

#!/bin/bash -l
#SBATCH --job-name=<jobname>
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=<sunetid>@stanford.edu
#SBATCH --output=<jobname>.%j.out
#SBATCH --error=<jobname>.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

## 3. Enter the lines of code to be run after the previous line

## 4. Press "esc" to stop editing. Type ":wq" to write and save. If required, you could also type ":q!" to quit without saving.

## 5. Run the sbatch file using the sbatch command

sbatch <sbatch_filename>.sbatch

## Note: For running sbatch files, the maximum time you can request is 48 hours. This might result in the job/process getting timed out before it is completed. Split your tasks accordingly and submit as smaller jobs if required. For about 10 - 50 Gb Gb of fast5's, I typically run guppy_basecaller and guppy_barcoder within a few hours. Megalodon takes a few hours - a day, depending on the number of barcodes. Can submit Megalodon jobs for individual barcodes (or subset of barcodes), instead of a for loop over all the barcodes to avoid timing out.

## The following are the commands/scripts for the workflow:

# For fast basecalling - fast5 to fastq

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '<input_directory_containing_fast5s>' --save_path '<output_directory_for_fastq>' --config '/home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg' --num_callers 32 --fast5_out

# Splitting reads in fastq based on barcode - "demultiplexing"

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_barcoder --input_path '<output_directory_for_fastq>/pass' --save_path '<output_directory_for_fastq>/split' --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/barcoding/configuration.cfg -t 30

# Reading readid's from each fastq file and writing into a txt file "barcode01_readids.txt" and so on - for each barcode - this txt file will be used by "fast5_subset" to split/demultiplex the fast5 file by barcodes

cd  <output_directory_for_fastq>/split

#the following is for barcodes 01 - 24, 
for bar in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
do
	for file in './barcode'$bar'/*fastq';
	do
	awk -F '[ ]' 'BEGIN {OFS="\n"}(NR%4==1){print $1}' $file | sed -e 's/^.//' >> './barcode'$bar'/barcode'$bar'_readids.txt'
	done
done

#number of reads per barcode can be output using the following line of code (wc -l <filename> reads the number of lines in the file)

wc -l barcode*/*.txt

#Splitting fast5s by barcode
#Making a new folder fast5_split

cd <output_directory_for_fastq>

mkdir fast5_split

#Making a new folder for each barcode within fast5_split

for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
do
mkdir 'fast5_split/barcode'$n
done

#activate the environment that contains fast5_subset. Running fast5_subset which reads a list of readid's and output the corresponding reads into smaller/subset fast5 files

conda activate ont2

for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
do
fast5_subset -i '<input_directory_containing_fast5s>' -s '<output_directory_for_fastq>/fast5_split/barcode'$n -l '<output_directory_for_fastq>/split/barcode'$n'/barcode'$n'_readids.txt' -n 2000000 -t 30 
done

###

#Running Megalodon for all barcodes

###
mkdir <megalodon_output_directory>

#Activate environment containing Megalodon
#activate samtools for modifying/reading bam files output by megalodon

conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools

#Running megalodon for all barcodes 01 - 24. This can also be done for each barcode individually
#<template_file>.fa refers to the fasta format file of the plasmid or reference template. For example, /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa


for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
do
megalodon '<output_directory_for_fastq>/fast5_split/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '<megalodon_output_directory>/barcode'$n --reference <template_file>.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' 'EXP-NBD114' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0
done

#It is recommended to filter files based on length close in range to the template length. This avoids errors in reading the bam files later in Jupyter. This will get rid of outliers that are much longer or shorter than template etc. I typically do +/-~100 bp to the expected read length. For example 450 - 600, for a ~500 bp expected read length.

for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
do
/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='<megalodon_output_directory>/barcode'$n'/mod_mappings.bam' out='<megalodon_output_directory>/barcode'$n'/mod_mappings_<approx_length>bp.bam' minlength=<min_length_in_bp> maxlength=<max_length_in_bp> ref=<template_file>.fa
done

######
## Megalodon generated mod_mappings_<approx_length>bp.bam files are the required output for running Jupyter notebooks to read mA (and mC) and for dimelo data analyses.
## The bam files can be viewed with samtools to make sure there are enough reads in each barcode as follows

conda activate ont2
ml biology samtools

#for all reads output by megalodon
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
do
samtools view '<megalodon_output_directory>/barcode'$n'/mod_mappings.bam' | wc -l
done


#for length filtered reads output by megalodon and then bbmap's reformat.sh
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24;
do
samtools view '<megalodon_output_directory>/barcode'$n'/mod_mappings_<approx_length>bp.bam' | wc -l
done
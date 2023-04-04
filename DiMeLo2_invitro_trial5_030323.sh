#DiMeLo2_invitro_trial5
#02.24.23
#Workflow:
#Fast basecall using Guppy fast mode
#Split barcodes using fastq
#separate fast5 based on split fastq
#megalodon mod basecall

#####
ssh sh02-ln04
screen -S invitro

#####

#Making backup:

#mkdir /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro
cd /oak/stanford/groups/astraigh/minION/all_data_backup/data

tar --use-compress-program="pigz" -cvf DiMeLo_v2/invitro/DiMeLo2__invitro_trial5_2023-03-05.tar.gz DiMeLo_v2_invitro_trial5 >> DiMeLo_v2/invitro/tarlog_2023-03-05.txt

ml system rclone

rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro/DiMeLo2_invitro_trial5_2023-03-05.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/invitro" -v --dry-run

#####

#For fast basecalling (we use this mainly to demultiplex)

#Did this part also on GPU node

#sbatch file location: 
vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial5_fast.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial5_bar_fast.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial5_bar_fast.guppy.full.%j.out
#SBATCH --error=DiMeLo_v2_invitro_trial5_bar_fast.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial5/v2_trial5/20230303_1222_MN28439_FAV93080_d36b55af/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial5/full_fast_gpu' --config '/home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg' --num_callers 32 --fast5_out

#submitted with the ^ sbatch script (hopefully runs within 48 hours total time)

#demux using guppy_barcoder
#salloc -p gpu -c 30 -G 2 --time=48:00:00
#ssh sh03-12n14 #or w.e node you get

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_barcoder --input_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial5/full_fast_gpu/pass --save_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial5/full_fast_gpu/split --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/barcoding/configuration.cfg -t 30

#making readid list by barcode

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial5/full_fast_gpu/split

for bar in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
	for file in './barcode'$bar'/*fastq';
	do
	awk -F '[ ]' 'BEGIN {OFS="\n"}(NR%4==1){print $1}' $file | sed -e 's/^.//' >> './barcode'$bar'/barcode'$bar'_readids.txt'
	done
done

wc -l barcode*/*.txt

#   5406 barcode01/barcode01_readids.txt
#   6835 barcode02/barcode02_readids.txt
#   6287 barcode03/barcode03_readids.txt
#  12960 barcode04/barcode04_readids.txt
#   4491 barcode05/barcode05_readids.txt
#   5761 barcode06/barcode06_readids.txt
#  21731 barcode07/barcode07_readids.txt
#   1084 barcode08/barcode08_readids.txt
#   8446 barcode09/barcode09_readids.txt
#  13064 barcode10/barcode10_readids.txt
#  11368 barcode11/barcode11_readids.txt
#   5880 barcode12/barcode12_readids.txt
#   9501 barcode13/barcode13_readids.txt
#   5674 barcode14/barcode14_readids.txt
#   3782 barcode15/barcode15_readids.txt
#  10406 barcode16/barcode16_readids.txt
# 132676 total

#Splitting fast5s by barcode
#in gpu node

conda activate ont2
cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial5/full_fast_gpu

mkdir fast5_split

for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
mkdir 'fast5_split/barcode'$n
done
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
fast5_subset -i '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial5/v2_trial5/20230303_1222_MN28439_FAV93080_d36b55af/fast5' -s '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial5/full_fast_gpu/fast5_split/barcode'$n -l '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial5/full_fast_gpu/split/barcode'$n'/barcode'$n'_readids.txt' -n 2000000 -t 30 
done

###

 10997 barcode01/filename_mapping.txt
 31242 barcode02/filename_mapping.txt
 20865 barcode03/filename_mapping.txt
 32827 barcode04/filename_mapping.txt
  9122 barcode05/filename_mapping.txt
 13793 barcode06/filename_mapping.txt
 24708 barcode07/filename_mapping.txt
 40610 barcode08/filename_mapping.txt
 38819 barcode09/filename_mapping.txt
 28497 barcode10/filename_mapping.txt
 40451 barcode11/filename_mapping.txt
 19969 barcode12/filename_mapping.txt
 32176 barcode13/filename_mapping.txt
 14041 barcode14/filename_mapping.txt
 36380 barcode15/filename_mapping.txt
 24042 barcode16/filename_mapping.txt
 23651 barcode17/filename_mapping.txt
 48354 barcode18/filename_mapping.txt
490544 total


#Running Megalodon for all barcodes
#Running everything in GPU node (as number of reads is small)


###
mkdir /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial5/

vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial5.megalodon.samtools.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial5_bar.megalodon.samtools
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial5_bar.megalodon.samtools.%j.out
#SBATCH --error=DiMeLo_v2_invitro_trial5_bar.megalodon.samtools.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g



conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
megalodon '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial5/full_fast_gpu/fast5_split/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial5/barcode'$n --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial5/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial5/barcode'$n'/mod_mappings_500bp.bam' minlength=450 maxlength=600 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial5/barcode'$n'/mod_mappings.bam' | wc -l >> '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial5/readcounts.txt'
done

for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial5/barcode'$n'/mod_mappings_500bp.bam' | wc -l '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial5/readcounts_500bp.txt'
done


#DiMeLo2_invitro_trial2
#11.10.22
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

mkdir /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro
cd /oak/stanford/groups/astraigh/minION/all_data_backup/data

tar --use-compress-program="pigz" -cvf DiMeLo_v2/invitro/DiMeLo2__invitro_trial2_2022-11-10.tar.gz DiMeLo_v2_invitro_trial2 >> DiMeLo_v2/invitro/tarlog_2022-11-10.txt

ml system rclone
rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro/DiMeLo2__invitro_trial2_2022-11-10.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/invitro" -v --dry-run

rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro/DiMeLo2__invitro_trial2_2022-11-10.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/invitro" -v

#####

#For fast basecalling (we use this mainly to demultiplex)

#sbatch file location: 
vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial2_fast.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial2_bar_fast.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial2_bar_fast.guppy.full.%j.out
#SBATCH --error=DiMeLo_v2_invitro_trial2_bar_fast.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial2/v2_trial2/20221109_1341_MN28439_FAU43818_2755e119/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial2/full_fast_gpu' --config '/home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg' --num_callers 32 --fast5_out

#submitted with the ^ sbatch script (hopefully runs within 48 hours total time)

#demux using guppy_barcoder
#salloc -p gpu -c 30 -G 2 --time=48:00:00
#ssh sh03-12n14 #or w.e node you get

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_barcoder --input_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial2/full_fast_gpu/pass --save_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial2/full_fast_gpu/split --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/barcoding/configuration.cfg -t 30

#making readid list by barcode

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial2/full_fast_gpu/split

for bar in 21 22 03 04 23 06 07 08 09 10;
do
	for file in './barcode'$bar'/*fastq';
	do
	awk -F '[ ]' 'BEGIN {OFS="\n"}(NR%4==1){print $1}' $file | sed -e 's/^.//' >> './barcode'$bar'/barcode'$bar'_readids.txt'
	done
done

wc -l barcode*/*.txt

###
    7881 barcode21/barcode21_readids.txt
    22308 barcode22/barcode22_readids.txt
    14410 barcode03/barcode03_readids.txt
    17290 barcode04/barcode04_readids.txt
    33634 barcode23/barcode23_readids.txt
    31797 barcode06/barcode06_readids.txt
    43096 barcode07/barcode07_readids.txt
    5869 barcode08/barcode08_readids.txt
    16211 barcode09/barcode09_readids.txt
    16955 barcode10/barcode10_readids.txt

   209451 total

###

#Splitting fast5s by barcode
#in gpu node

conda activate ont2
cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial2/full_fast_gpu

mkdir fast5_split

for n in 21 22 03 04 23 06 07 08 09 10;
do
mkdir 'fast5_split/barcode'$n
done
for n in 21 22 03 04 23 06 07 08 09 10;
do
fast5_subset -i '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial2/v2_trial2/20221109_1341_MN28439_FAU43818_2755e119/fast5' -s '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial2/full_fast_gpu/fast5_split/barcode'$n -l '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial2/full_fast_gpu/split/barcode'$n'/barcode'$n'_readids.txt' -n 2000000 -t 30 
done

###

#Running Megalodon for all barcodes
#Running everything in GPU node (as number of reads is small)

mkdir /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/

###

conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools
mkdir /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/
for n in 21 22 03 04 23 06 07 08 09 10;
do
megalodon '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial2/full_fast_gpu/fast5_split/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings_730bp.bam' minlength=700 maxlength=850 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 21 22 03 04 23 06 07 08 09 10;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings.bam' | wc -l
done

for n in 21 22 03 04 23 06 07 08 09 10;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings_730bp.bam' | wc -l
done

'''
mod_mappings.bam    mod_mappings_730bp.bam
5718                3153    
15864               5243
7745                4527
21931               7697
16600               11858
13821               10457
17011               12770
32985               26676
31338               24075
42042               33634
'''



###Same thing as above but in sbatch

vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial2_megalodon.samtools.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial2_bar.megalodon.samtools
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial2_bar.megalodon.samtools.%j.out
#SBATCH --error=DiMeLo_v2_invitro_trial2_bar.megalodon.samtools.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools
#mkdir /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/
for n in 21 22 03 04 23 06 07 08 09 10;
do
megalodon '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial2/full_fast_gpu/fast5_split/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n --overwrite --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings_730bp.bam' minlength=700 maxlength=850 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 21 22 03 04 23 06 07 08 09 10;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings.bam' | wc -l
done

for n in 21 22 03 04 23 06 07 08 09 10;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings_730bp.bam' | wc -l
done

#Some reason no bam file in barcode09, repeating megalodon for this folder

vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial2_megalodon.samtools.09.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial2_bar.megalodon.samtools.09
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial2_bar.megalodon.samtools.%j.09.out
#SBATCH --error=DiMeLo_v2_invitro_trial2_bar.megalodon.samtools.%j.09.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools
#mkdir /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/
for n in 09;
do
megalodon '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial2/full_fast_gpu/fast5_split/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n --overwrite --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings_730bp.bam' minlength=700 maxlength=850 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 09;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings.bam' | wc -l
done

for n in 09;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial2/barcode'$n'/mod_mappings_730bp.bam' | wc -l
done


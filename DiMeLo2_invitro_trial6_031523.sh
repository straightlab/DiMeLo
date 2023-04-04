#DiMeLo2_invitro_trial6
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

tar --use-compress-program="pigz" -cvf DiMeLo_v2/invitro/DiMeLo2_invitro_trial6_2023-03-29.tar.gz DiMeLo_v2_invitro_trial6 >> DiMeLo_v2/invitro/tarlog_2023-03-29.txt

ml system rclone

rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro/DiMeLo2__invitro_trial6_2023-03-29.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/invitro" -v

#####

#For fast basecalling (we use this mainly to demultiplex)

#Did this part also on GPU node

#sbatch file location: 
vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial6_fast.guppy_demux.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial6_bar_fast.guppy_demux.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial6_bar_fast.guppy_demux.full.%j.out
#SBATCH --error=DiMeLo_v2_invitro_trial6_bar_fast.guppy_demux.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial6/v2_trial6/20230313_1215_MN28439_FAV96273_60a63824/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu' --config '/home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg' --num_callers 32 --fast5_out


/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_barcoder --input_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu/pass --save_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu/split --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/barcoding/configuration.cfg -t 30

#making readid list by barcode

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu/split

for bar in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
	for file in './barcode'$bar'/*fastq';
	do
	awk -F '[ ]' 'BEGIN {OFS="\n"}(NR%4==1){print $1}' $file | sed -e 's/^.//' >> './barcode'$bar'/barcode'$bar'_readids.txt'
	done
done

wc -l barcode*/*.txt


#Splitting fast5s by barcode

conda activate ont2
cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu

mkdir fast5_split

for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
mkdir 'fast5_split/barcode'$n
done
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
fast5_subset -i '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial6/v2_trial6/20230313_1215_MN28439_FAV96273_60a63824/fast5' -s '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu/fast5_split/barcode'$n -l '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu/split/barcode'$n'/barcode'$n'_readids.txt' -n 2000000 -t 30 
done

###

###
#fast5_subset takes a while to run, so instead trying to run it on a smaller subset of all the reads that it is looking for in each barcode
#setting max number of reads per barcode to 10000

#sbatch file location: 
vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial6_fast.demux_10k.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial6_bar_fast.demux_10k.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial6_bar_fast.demux_10k.full.%j.out
#SBATCH --error=DiMeLo_v2_invitro_trial6_bar_fast.demux_10k.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

conda activate ont2
cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu
mkdir fast5_split_subset10k
cd fast5_split_subset10k
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
mkdir 'barcode'$n
done
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
head -10000 '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu/split/barcode'$n'/barcode'$n'_readids.txt' > '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu/split/barcode'$n'/barcode'$n'_readids_10k.txt'
done
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
fast5_subset -i '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial6/v2_trial6/20230313_1215_MN28439_FAV96273_60a63824/fast5' -s '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu/fast5_split_subset10k/barcode'$n -l '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu/split/barcode'$n'/barcode'$n'_readids_10k.txt' -n 2000000 -t 30 
done



#Running Megalodon for all barcodes
#Running everything in GPU node (as number of reads is small)


###
mkdir /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial6_10ksubset

vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial6.megalodon.samtools.10ksubset.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial6_bar.megalodon.samtools.10ksubset
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial6_bar.megalodon.samtools.10ksubset.%j.out
#SBATCH --error=DiMeLo_v2_invitro_trial6_bar.megalodon.samtools.10ksubset.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g



conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
megalodon '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial6/full_fast_gpu/fast5_split_subset10k/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial6_10ksubset/barcode'$n --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' 'EXP-NBD114' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial6_10ksubset/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial6_10ksubset/barcode'$n'/mod_mappings_500bp.bam' minlength=450 maxlength=600 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial6_10ksubset/barcode'$n'/mod_mappings.bam' | wc -l >> '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial6_10ksubset/readcounts.txt'
done

for n in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial6_10ksubset/barcode'$n'/mod_mappings_500bp.bam' | wc -l >> '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial6_10ksubset/readcounts_500bp.txt'
done


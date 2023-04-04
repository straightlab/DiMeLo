#DiMeLo2_invitro_trial3
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

#mkdir /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro
cd /oak/stanford/groups/astraigh/minION/all_data_backup/data

tar --use-compress-program="pigz" -cvf DiMeLo_v2/invitro/DiMeLo2__invitro_trial3_2022-11-17.tar.gz DiMeLo_v2_invitro_trial3 >> DiMeLo_v2/invitro/tarlog_2022-11-17.txt

ml system rclone

rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro/DiMeLo2__invitro_trial3_2022-11-17.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/invitro" -v

#####

#For fast basecalling (we use this mainly to demultiplex)

#sbatch file location: 
vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial3_fast.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial3_bar_fast.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial3_bar_fast.guppy.full.%j.out
#SBATCH --error=DiMeLo_v2_invitro_trial3_bar_fast.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial3/v2_trial3/20221117_1413_MN28439_FAU63227_444609f3/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu' --config '/home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg' --num_callers 32 --fast5_out

#submitted with the ^ sbatch script (hopefully runs within 48 hours total time)

#demux using guppy_barcoder
#salloc -p gpu -c 30 -G 2 --time=48:00:00
#ssh sh03-12n14 #or w.e node you get

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_barcoder --input_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu/pass --save_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu/split --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/barcoding/configuration.cfg -t 30

#making readid list by barcode

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu/split

for bar in 21 22 03 04 23 24 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
	for file in './barcode'$bar'/*fastq';
	do
	awk -F '[ ]' 'BEGIN {OFS="\n"}(NR%4==1){print $1}' $file | sed -e 's/^.//' >> './barcode'$bar'/barcode'$bar'_readids.txt'
	done
done

wc -l barcode*/*.txt

###

#21    8834
#22   25017
#03   94400
#04   12242
#23   22261
#24  148539
#07   14152
#08   19274
#09  148937
#10   33699
#11   57430
#12   20158
#13  186626
#14   22995
#15  178672
#16   24884
#17  151513
#18  243688
#19  165057
#20  113679
#----------
#sum 209451

###

#Splitting fast5s by barcode
#in gpu node

conda activate ont2
cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu

mkdir fast5_split

for n in 21 22 03 04 23 24 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
mkdir 'fast5_split/barcode'$n
done
for n in 21 22 03 04 23 24 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
fast5_subset -i '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial3/v2_trial3/20221117_1413_MN28439_FAU63227_444609f3/fast5' -s '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu/fast5_split/barcode'$n -l '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu/split/barcode'$n'/barcode'$n'_readids.txt' -n 2000000 -t 30 
done

#Got interrupted at 96% for 16, ran 17 - 20 in new job.
###

#Running Megalodon for all barcodes
#Running everything in GPU node (as number of reads is small)


###
mkdir /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/

conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools
for n in 21 22 03 04 23 24 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
megalodon '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu/fast5_split/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings_730bp.bam' minlength=700 maxlength=850 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 21 22 03 04 23 24 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings.bam' | wc -l
done

for n in 21 22 03 04 23 24 07 08 09 10 11 12 13 14 15 16 17 18 19 20;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings_730bp.bam' | wc -l
done



###Megalodon onwards but in sbatch



vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial3.megalodon.samtools.1.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.%j.1.out
#SBATCH --error=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.%j.1.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools
for n in 21 22 03 04 23;
do
megalodon '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu/fast5_split/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n --overwrite --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings_730bp.bam' minlength=700 maxlength=850 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 21 22 03 04 23;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings.bam' | wc -l
done

for n in 21 22 03 04 23;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings_730bp.bam' | wc -l
done

vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial3.megalodon.samtools.2.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.2
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.%j.2.out
#SBATCH --error=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.%j.2.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools
for n in 24 07 08 09 10;
do
megalodon '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu/fast5_split/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n --overwrite --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings_730bp.bam' minlength=700 maxlength=850 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 24 07 08 09 10;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings.bam' | wc -l
done

for n in 24 07 08 09 10;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings_730bp.bam' | wc -l
done

vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial3.megalodon.samtools.3.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.3
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.%j.3.out
#SBATCH --error=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.%j.3.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools
for n in 11 12 13 14 15;
do
megalodon '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu/fast5_split/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n --overwrite --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings_730bp.bam' minlength=700 maxlength=850 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 11 12 13 14 15;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings.bam' | wc -l
done

for n in 11 12 13 14 15;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings_730bp.bam' | wc -l
done

vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial3.megalodon.samtools.4.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.%j.4.out
#SBATCH --error=DiMeLo_v2_invitro_trial3_bar.megalodon.samtools.%j.4.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

conda activate /home/groups/astraigh/miniconda3/envs/ont/
ml devel java
ml biology samtools
for n in 16 17 18 19 20;
do
megalodon '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial3/full_fast_gpu/fast5_split/barcode'$n --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n --overwrite --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings_730bp.bam' minlength=700 maxlength=850 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 16 17 18 19 20;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings.bam' | wc -l
done

for n in 16 17 18 19 20;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial3/barcode'$n'/mod_mappings_730bp.bam' | wc -l
done
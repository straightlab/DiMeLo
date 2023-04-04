#DiMeLo2_invitro_trial1
#10.18.22
#Workflow:
#Fast basecall using Guppy fast mode
#Split barcodes using fastq
#separate fast5 based on split fastq
#megalodon mod basecall

#####
ssh sh03-ln08
screen -S invitro

#####

#Making backup:

mkdir /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro
cd /oak/stanford/groups/astraigh/minION/all_data_backup/data

tar --use-compress-program="pigz" -cvf DiMeLo_v2/invitro/DiMeLo2__invitro_trial1_2022-10-18.tar.gz DiMeLo_v2_invitro_trial1 >> DiMeLo_v2/invitro/tarlog_2022-10-18.txt

ml system rclone
rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro/DiMeLo2__invitro_trial1_2022-10-18.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/invitro" -v --dry-run

rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/invitro/DiMeLo2__invitro_trial1_2022-10-18.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/invitro" -v

#####

#For fast basecalling (we use this mainly to demultiplex)

#sbatch file location: 
vim ~/sbatch_scripts/DiMeLo_v2_invitro_trial1_fast.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo_v2_invitro_trial1_bar_fast.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo_v2_invitro_trial1_bar_fast.guppy.full.%j.out
#SBATCH --error=DiMeLo_v2_invitro_trial1_bar_fast.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial1/v2_trial1/20221013_1336_MN28439_FAU43818_0c2eb966/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu' --config '/home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg' --num_callers 32 --fast5_out

#10.20.22

#demux using guppy_barcoder
#salloc -p gpu -c 30 -G 2 --time=48:00:00
#ssh sh03-12n14 #or w.e node you get

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_barcoder --input_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu/pass --save_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu/split --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/barcoding/configuration.cfg -t 30

#making readid list by barcode

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu/split

for bar in 21 22 03 04 23 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
	for file in './barcode'$bar'/*fastq';
	do
	awk -F '[ ]' 'BEGIN {OFS="\n"}(NR%4==1){print $1}' $file | sed -e 's/^.//' >> './barcode'$bar'/barcode'$bar'_readids.txt'
	done
done

wc -l barcode*/*.txt

###
    11741 barcode21/barcode21_readids.txt
   101128 barcode22/barcode22_readids.txt
    12333 barcode03/barcode03_readids.txt
    89136 barcode04/barcode04_readids.txt
    10051 barcode23/barcode23_readids.txt
    72519 barcode06/barcode06_readids.txt
     7817 barcode07/barcode07_readids.txt
    14918 barcode08/barcode08_readids.txt
    33917 barcode09/barcode09_readids.txt
    47863 barcode10/barcode10_readids.txt
    51848 barcode11/barcode11_readids.txt
    47731 barcode12/barcode12_readids.txt
    36986 barcode13/barcode13_readids.txt
    49359 barcode14/barcode14_readids.txt
    29170 barcode15/barcode15_readids.txt
    49048 barcode16/barcode16_readids.txt
    62468 barcode17/barcode17_readids.txt
    70213 barcode18/barcode18_readids.txt
   798246 total

###

#Splitting fast5s by barcode

conda activate ont2
cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu

mkdir fast5_split

for n in 21 22 03 04 23 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
mkdir 'fast5_split/barcode'$n
done
for n in 21 22 03 04 23 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
fast5_subset -i '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2_invitro_trial1/v2_trial1/20221013_1336_MN28439_FAU43818_0c2eb966/fast5' -s '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu/fast5_split/barcode'$n -l '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu/split/barcode'$n'/barcode'$n'_readids.txt' -n 2000000 -t 30 
done

###

#Running Megalodon for all barcodes
#barcode 11 as template

mkdir /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/

vim /home/groups/astraigh/ont_minion/scripts/sbatch/D2_invitro_trial1_1x601.megalodon_barcode11_gpu.sbatch

---

#!/bin/bash -l
#SBATCH --job-name=D2_invitro_trial1_1x601.megalodon.barcode11.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=D2_invitro_trial1_1x601.megalodon.barcode11.full.%j.out
#SBATCH --error=D2_invitro_trial1_1x601.megalodon.barcode11.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
conda activate /home/groups/astraigh/miniconda3/envs/ont/
megalodon /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu/fast5_split/barcode11 --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/barcode11 --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

megalodon /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu/fast5_split/barcode11 --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/barcode11_mC --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases_5mC_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

---

sbatch /home/groups/astraigh/ont_minion/scripts/sbatch/D2_invitro_trial1_1x601.megalodon_barcode11_gpu.sbatch

###
#Using for loop to make megalodon sbatch scripts for all barcodes


for n in 21 22 03 04 23 06 07 08 09 10 12 13 14 15 16 17 18;
do
sed 's/barcode11/barcode'$n'/g' D2_invitro_trial1_1x601.megalodon_barcode11_gpu.sbatch > 'D2_invitro_trial1_1x601.megalodon_barcode'$n'_gpu.sbatch'
sbatch 'D2_invitro_trial1_1x601.megalodon_barcode'$n'_gpu.sbatch'
done

###


ml devel java
cd /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/
for n in 21 22 03 04 23 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/barcode'$n'/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/barcode'$n'/mod_mappings_730bp.bam' minlength=700 maxlength=850 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

ml devel java
cd /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/
for n in 21 22 03 04 23 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
/home/groups/astraigh/miniconda3/envs/charseq/opt/bbmap-38.06/reformat.sh in='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/barcode'$n'_mC/mod_mappings.bam' out='/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/barcode'$n'_mC/mod_mappings_730bp.bam' minlength=700 maxlength=850 ref=/home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa
done

conda activate ont2
ml biology samtools
for n in 21 22 03 04 23 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/barcode'$n'/mod_mappings.bam' | wc -l
done

for n in 21 22 03 04 23 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/barcode'$n'/mod_mappings_730bp.bam' | wc -l
done

'''
7928
79388
8498
69581
6450
57083
4802
9623
17839
15005
37640
35832
27885
0
20743
0
0
0
'''

for n in 21 22 03 04 23 06 07 08 09 10 11 12 13 14 15 16 17 18;
do
samtools view '/scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/barcode'$n'_mC/mod_mappings_730bp.bam' | wc -l
done

'''
4228
75741
4599
68320
3560
56293
2344
5847
12632
11212
37162
35160
27484
36147
20575
35751
33649
21481
'''


vim /home/groups/astraigh/ont_minion/scripts/sbatch/D2_invitro_trial1_1x601.megalodon_barcode14_gpu_redone.sbatch

#!/bin/bash -l
#SBATCH --job-name=D2_invitro_trial1_1x601.megalodon.barcode14_redone.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=D2_invitro_trial1_1x601.megalodon.barcode14_redone.full.%j.out
#SBATCH --error=D2_invitro_trial1_1x601.megalodon.barcode14_redone.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
conda activate /home/groups/astraigh/miniconda3/envs/ont/
megalodon /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu/fast5_split/barcode14 --outputs basecalls mod_basecalls mappings mod_mappings mods per_read_mods --output-directory /scratch/groups/astraigh/minion_seq/megalodon/DiMeLo_v2_invitro_trial1/barcode14_redone --reference /home/groups/astraigh/ont_minion/plasmids/pGEM_3z_601.fa --guppy-server-path /home/groups/astraigh/software/guppy_4.4.2_gpu/ont-guppy/bin/guppy_basecall_server --guppy-config res_dna_r941_min_modbases-all-context_v001.cfg --guppy-params "-d /home/groups/astraigh/software/rerio/basecall_models/ --num_callers 30 --barcode_kits 'EXP-NBD104' --trim_barcodes" --guppy-timeout 360 --mod-min-prob 0

sbatch /home/groups/astraigh/ont_minion/scripts/sbatch/D2_invitro_trial1_1x601.megalodon_barcode14_gpu_redone.sbatch

mkdir /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/guppy_allcontext

###All context guppy modbasecalling without alignment
vim /home/groups/astraigh/ont_minion/scripts/sbatch/D2_invitro_trial1_1x601.guppy_allcontext_barcode14_gpu.sbatch

#!/bin/bash -l
#SBATCH --job-name=D2_invitro_trial1_1x601.guppy_allcontext_barcode14_gpu.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=D2_invitro_trial1_1x601.guppy_allcontext_barcode14_gpu.full.%j.out
#SBATCH --error=D2_invitro_trial1_1x601.guppy_allcontext_barcode14_gpu.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu/fast5_split/barcode14/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/guppy_allcontext/barcode14' -c res_dna_r941_min_modbases-all-context_v001.cfg -d /home/groups/astraigh/software/rerio/basecall_models/ --bam_out --bam_methylation_threshold 0.0 --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 4096 --max_queued_reads 6000 -x 'auto'

sbatch /home/groups/astraigh/ont_minion/scripts/sbatch/D2_invitro_trial1_1x601.guppy_allcontext_barcode14_gpu.sbatch


cd /home/groups/astraigh/ont_minion/scripts/sbatch/
for n in 23 11;
do
sed 's/barcode14/barcode'$n'/g' D2_invitro_trial1_1x601.guppy_allcontext_barcode14_gpu.sbatch > 'D2_invitro_trial1_1x601.guppy_allcontext_barcode'$n'_gpu.sbatch'
sbatch 'D2_invitro_trial1_1x601.guppy_allcontext_barcode'$n'_gpu.sbatch'
done



mkdir /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/guppy_mC

###mC guppy modbasecalling without alignment
vim /home/groups/astraigh/ont_minion/scripts/sbatch/D2_invitro_trial1_1x601.guppy_mC_barcode14_gpu_lowload.sbatch

#!/bin/bash -l
#SBATCH --job-name=D2_invitro_trial1_1x601.guppy_mC_barcode14_gpu.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=D2_invitro_trial1_1x601.guppy_mC_barcode14_gpu.full.%j.out
#SBATCH --error=D2_invitro_trial1_1x601.guppy_mC_barcode14_gpu.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/full_fast_gpu/fast5_split/barcode14/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v2_invitro_trial1/guppy_mC/barcode14' -c res_dna_r941_min_modbases_5mC_v001.cfg -d /home/groups/astraigh/software/rerio/basecall_models/ --bam_out --bam_methylation_threshold 0.0 --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 2048 --max_queued_reads 3000 -x 'auto'


sbatch /home/groups/astraigh/ont_minion/scripts/sbatch/D2_invitro_trial1_1x601.guppy_mC_barcode14_gpu_lowload.sbatch


cd /home/groups/astraigh/ont_minion/scripts/sbatch/
for n in 23 11;
do
sed 's/barcode14/barcode'$n'/g' D2_invitro_trial1_1x601.guppy_mC_barcode14_gpu_lowload.sbatch > 'D2_invitro_trial1_1x601.guppy_mC_barcode'$n'_gpu_lowload.sbatch'
sbatch 'D2_invitro_trial1_1x601.guppy_mC_barcode'$n'_gpu_lowload.sbatch'
done
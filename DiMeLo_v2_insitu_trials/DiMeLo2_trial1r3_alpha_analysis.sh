#7/04/22

#Making backups first
mkdir /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2
mkdir /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/insitu
cd /oak/stanford/groups/astraigh/minION/all_data_backup/data

tar --use-compress-program="pigz" -cvf DiMeLo_v2/insitu/DiMeLo2_trial1r3_alpha_2022-07-04.tar.gz DiMeLo2_trial1r3_alpha >> DiMeLo_v2/insitu/tarlog_2022-06-23.txt

ml system rclone
rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/insitu/DiMeLo2_trial1r3_alpha_2022-07-04.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/insitu" -v --dry-run

rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/insitu/DiMeLo2_trial1r3_alpha_2022-07-04.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/insitu" -v

######

#6/23/22

#For fast basecalling (we use this mainly to demultiplex)

#sbatch file location: ~/sbatch_scripts/DiMeLo2_trial1r3_alpha_fast.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1r3_alpha_bar_fast.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1r3_alpha_bar_fast.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1r3_alpha_bar_fast.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo2_trial1r3_alpha/D2_trial1r3_alpha/20220630_1449_MN28439_FAT21138_87fcde79/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu' --config '/home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg' --num_callers 32 --fast5_out


#06/25/22

#first demultiplex fastq using guppy_barcoder
#submit as a sbatch

#if not sbatch
#salloc -p gpu -c 30 -G 2 --time=48:00:00
#ssh sh03-12n09 #or w.e node you get


#sbatch file location: ~/sbatch_scripts/DiMeLo2_trial1r3_alpha_demux.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1r3_alpha_bar_fast_demux.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1r3_alpha_bar_fast.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1r3_alpha_bar_fast.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_barcoder --input_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu/pass --save_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu/split --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/barcoding/configuration.cfg -t 30

sbatch ~/sbatch_scripts/DiMeLo2_trial1r3_alpha_demux.guppy.sbatch

---

#making a list of all readids for each barcode

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu/split

for bar in 01 02 05 06;
do
	for file in './barcode'$bar'/*fastq';
	do
	awk -F '[ ]' 'BEGIN {OFS="\n"}(NR%4==1){print $1}' $file | sed -e 's/^.//' >> './barcode'$bar'/barcode'$bar'_readids.txt'
	done
done

'''
#optional
wc -l /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu/split/barcode*/*.txt
291475 barcode01/barcode01_readids.txt
586081 barcode02/barcode02_readids.txt
549753 barcode05/barcode05_readids.txt
570911 barcode06/barcode06_readids.txt
1998220 total 

'''

######

##split raw fast5s by barcode readids to run through megalodon
conda activate ont2

###need to mkdir for outputdirs ./fast5_split and for each ./fast5_split/barcode0*
cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu

mkdir fast5_split

for n in 01 02 05 06;
do
mkdir 'fast5_split/barcode'$n
done


for n in 01 02 05 06;
do
fast5_subset -i /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo2_trial1r3_alpha/D2_trial1r3_alpha/20220630_1449_MN28439_FAT21138_87fcde79/fast5 -s '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu/fast5_split/barcode'$n -l '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu/split/barcode'$n'/barcode'$n'_readids.txt' -n 4000 -t 32 
done


'''
#optional
wc -l /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu/fast5_split/*/filename_mapping.txt

291475 barcode01/filename_mapping.txt
586081 barcode02/filename_mapping.txt
549753 barcode05/filename_mapping.txt
570911 barcode06/filename_mapping.txt
1998220 total

'''

######

###
#superaccuracy basecalling for better quality mapping

#sbatch file location: ~/sbatch_scripts/DiMeLo2_trial1r3_alpha_bar_supacc.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1r3_alpha_bar_supacc.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1r3_alpha_bar_supacc.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1r3_alpha_bar_supacc.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo2_trial1r3_alpha/D2_trial1r3_alpha/20220630_1449_MN28439_FAT21138_87fcde79/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha_sup/' -c dna_r9.4.1_450bps_sup.cfg -d /home/groups/astraigh/software/guppy_tar/ont-guppy/data/ --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 128 --max_queued_reads 6000 -x 'auto'


###
#for all context mod basecalling after splitting fast5s

#sbatch file location: ~/sbatch_scripts/DiMeLo2_trial1r3_alpha_bar_allcontext.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1r3_alpha_bar_allcontext.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1r3_alpha_bar_allcontext.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1r3_alpha_bar_allcontext.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo2_trial1r3_alpha/D2_trial1r3_alpha/20220630_1449_MN28439_FAT21138_87fcde79/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha_mod/' -c res_dna_r941_min_modbases-all-context_v001.cfg -d /home/groups/astraigh/software/rerio/basecall_models/ --bam_out --bam_methylation_threshold 0.0 --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 4096 --max_queued_reads 6000 -x 'auto'

######

#07/05/22

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha

for n in 01 02 05 06;
do
mkdir '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode'$n'_sup'
done

for n in 01 02 05 06;
do
mkdir '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode'$n'_mod'
done

#sbatch scripts or superacc basecalling and all_context mod basecalling of demultiplexed libraries:

#example:
#barcode06

#Supacc basecalling

#sbatch file location: 
vim ~/sbatch_scripts/DiMeLo2_trial1r3_alpha_barcode06_supacc.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1r3_alpha_barcode06_supacc.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1r3_alpha_barcode06_supacc.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1r3_alpha_barcode06_supacc.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu/fast5_split/barcode06/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode06_sup' -c dna_r9.4.1_450bps_sup.cfg -d /home/groups/astraigh/software/guppy_tar/ont-guppy/data/ --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 128 --max_queued_reads 6000 -x 'auto'

sbatch ~/sbatch_scripts/DiMeLo2_trial1r3_alpha_barcode06_supacc.guppy.sbatch

###
#for all context mod basecalling after splitting fast5s

#sbatch file location: 

vim ~/sbatch_scripts/DiMeLo2_trial1r3_alpha_barcode06_allcontext.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1r3_alpha_barcode06_allcontext.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1r3_alpha_barcode06_allcontext.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1r3_alpha_barcode06_allcontext.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/full_fast_gpu/fast5_split/barcode06/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode06_mod' -c res_dna_r941_min_modbases-all-context_v001.cfg -d /home/groups/astraigh/software/rerio/basecall_models/ --bam_out --bam_methylation_threshold 0.0 --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 4096 --max_queued_reads 6000 -x 'auto'

sbatch ~/sbatch_scripts/DiMeLo2_trial1r3_alpha_barcode06_allcontext.guppy.sbatch

#######

#merging modbasecalled fastq's in each barcode into 1 file

for bar in 01 02 05 06;
do
	cd '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode'$bar'_mod/pass'
	for file in fastq_runid*.fastq;
	do
	cat $file >> 'barcode'$bar'_merge.fastq'
	done
done

'''
(base) [ksundar@sh02-09n13 /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha] (job 56385736) $ wc -l barcode*_mod/pass/*merge.fastq
 	1117488 barcode01_mod/pass/barcode01_merge.fastq
 	2282940 barcode02_mod/pass/barcode02_merge.fastq
 	2109396 barcode05_mod/pass/barcode05_merge.fastq
 	2221164 barcode06_mod/pass/barcode06_merge.fastq
 	7730988 total

'''

#######

#Chm13 + HG002 ChrX ChrY genome

#This part was done earlier on 2022.02.07

#chm13_T2T autosomes + MT + ChrX and ChrY v2.7 from HG002
#Downloaded latest frankengenome from https://s3-us-west-2.amazonaws.com/human-pangenomics/working/T2T/HG002XY/v2.1/v2.7.fasta locally
#Uploaded frankengenome to /oak/stanford/groups/astraigh/T2T_HG002_XY_v2.7

#preprocess kmers for T2T_HG002_XY_v2.7
cd /home/groups/astraigh/ont_minion/meryl

/home/groups/astraigh/software/Winnowmap/bin/meryl count k=15 output /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_merylDB /oak/stanford/groups/astraigh/T2T_HG002_XY_v2.7
/home/groups/astraigh/software/Winnowmap/bin/meryl print greater-than distinct=0.9998 /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_merylDB > /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_repetitive_k15.txt

#######

#Running winnomap on superaccuracy basecalled fastq files to map to the Chm13 HG002 ChrX ChrY frankengenome (details above)

## merge fastq files within each sup folder

for bar in 01 02 05 06;
do
	cd '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode'$bar'_sup/pass'
	for file in fastq_runid*.fastq;
	do
	cat $file >> 'barcode'$bar'_merge.fastq'
	done
done

'''
(base) [ksundar@sh02-09n13 /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha] (job 56385736) $ wc -l barcode*_sup/pass/*merge.fastq
	1138604 barcode01_sup/pass/barcode01_merge.fastq
	2266044 barcode02_sup/pass/barcode02_merge.fastq
	2141164 barcode05_sup/pass/barcode05_merge.fastq
	2223444 barcode06_sup/pass/barcode06_merge.fastq
	7769256 total
'''

#winnowmap > bam; remove duplicates, unmapped etc.; sort and index

#If needed:
mkdir /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/


module load biology
module load samtools

for bar in 06;
do
### KS new code: run the following
/home/groups/astraigh/software/Winnowmap/bin/winnowmap -W /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_repetitive_k15.txt -ax map-ont /oak/stanford/groups/astraigh/T2T_HG002_XY_v2.7/v2.7_HG002_w_ChrY.fasta '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode'$bar'_sup/pass/barcode'$bar'_merge.fastq' | samtools view -b > '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.bam'
done


for bar in 06;
do
#first remove unmapped (4), secondary (256), and supplemental (2048) alignments, 4+256+2048=2308
samtools view -b -F 2308 '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.bam' > /'scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.clean.bam'
#then sort and index
samtools sort '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.clean.bam' > '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.clean.sorted.bam'
samtools index '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.clean.sorted.bam'
done
######

# of cleaned reads (unique mappers) by barcode
# 281007 - barcode01
# 559404 - barcode02
# 528855 - barcode05
# 549685 - barcode06

######
#Did the following clean up immediately following winnowmap:
'''
#got to clean up the bam to avoid duplicates
module load biology
module load samtools

#first remove unmapped (4), secondary (256), and supplemental (2048) alignments, 4+256+2048=2308
for bar in 01 02 05 06;
do
samtools view -b -@ 24 -F 2308 '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.sorted.bam' > /'scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.sorted.clean.bam'
done
'''
#062822
#Not required: but generated one anyway (personal access token for github, 6/28/22 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)


#mkdir /scratch/groups/astraigh/minion_seq/dm/
mkdir /scratch/groups/astraigh/minion_seq/dm/DiMeLo2_trial1r3_alpha

#Running dimelo functions

conda activate dimelo

python
import dimelo as dm
for bar in ['05']:
	dm.qc_report('/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'+bar+'_merge.clean.sorted.bam', 'barcode'+bar, '/scratch/groups/astraigh/minion_seq/dm/DiMeLo2_trial1r3_alpha')


dm.qc_report('/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode06_merge.clean.sorted.bam', 'barcode06', '/scratch/groups/astraigh/minion_seq/dm/DiMeLo2_trial1r3_alpha')

exit()

######

#Making merged bam files in barcode*_mod/pass/ to sort before merging with alignment bam_out

ml biology samtools
for bar in 01 02 05 06;
do
	cd '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode'$bar'_mod/pass'
	for file in *.bam;
	do
	samtools view $file >> 'barcode'$bar'_merge.txt'
	done
done

#####

##The following is OKS's code for merging the info from sup_acc + winnowmap bam and the mod_basecalled bam, modified for the current library

mkdir /scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/

for n in 01 02 05 06;
do
	mkdir '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n
	join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.9,1.10,1.11,2.12,2.13 <(samtools view -@ 8 '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$n'_merge.clean.sorted.bam' | sort -k1) <(sort -k1 '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode'$n'_mod/pass/barcode'$n'_merge.txt') > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.txt'
	#remove reads with empty Mm and Ml to avoid problems downstream
	grep -v 'Mm:Z:A+a;C+m;' '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.txt' > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.txt'
	#Copy header from winnowmap output bam into sam
	samtools view -H '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1r3_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$n'_merge.clean.sorted.bam' | cat > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.sam'
	#concatenate cleaned up merged txt file into sam file with header
	cat '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.txt' >> '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.sam'
	#convert sam to bam, filter by quality score of 10
	samtools view -@30 -b -q10 '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.sam' > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.bam'
	#sort and index bam
	samtools sort -@30 '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.bam' > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.bam'
	samtools index -@30 '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.bam'
	#generate bigwigs to look at coverage,
	/home/groups/astraigh/miniconda3/envs/deepcrap/bin/bamCoverage -p 30 -bs 100 -b '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.bam' -o '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.100bp.bw'
done

###########

dm.plot_browser(["/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode01/barcode01_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.bam", "/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode02/barcode02_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.bam", "/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode05/barcode05_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.bam", "/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1r3_alpha/barcode06/barcode06_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.bam"], ["CENPA", "H3K9me3", "CTCF", "Cohesin"], "chrX_hg002:57450000-57750000", "A+CG", "/scratch/groups/astraigh/minion_seq/dm/DiMeLo2_trial1r3_alpha", static=False)

#Had issues running this


conda activate charseq2

vim ~/CTCF_consensus_seq.fasta
'''
>CTCF_consensus_seq_01
CCGCGAGGAGGCAG
>CTCF_consensus_seq_02
CCGCGAGGTGGCAG
>CTCF_consensus_seq_03
CCGCGAGGGGGCAG
>CTCF_consensus_seq_04
CCGCGAGGCGGCAG
>CTCF_consensus_seq_05
CCGCGTGGAGGCAG
>CTCF_consensus_seq_06
CCGCGTGGTGGCAG
>CTCF_consensus_seq_07
CCGCGTGGGGGCAG
>CTCF_consensus_seq_08
CCGCGTGGCGGCAG
>CTCF_consensus_seq_09
CCGCGGGGAGGCAG
>CTCF_consensus_seq_10
CCGCGGGGTGGCAG
>CTCF_consensus_seq_11
CCGCGGGGGGGCAG
>CTCF_consensus_seq_12
CCGCGGGGCGGCAG
>CTCF_consensus_seq_13
CCGCGCGGAGGCAG
>CTCF_consensus_seq_14
CCGCGCGGTGGCAG
>CTCF_consensus_seq_15
CCGCGCGGGGGCAG
>CTCF_consensus_seq_16
CCGCGCGGCGGCAG
'''


bbduk.sh in='/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode05_sup/pass/barcode05_merge.fastq' -qin=33 ref='~/CTCF_consensus_seq.fasta' outm='/scratch/groups/astraigh/kousik/hg002_alignments/bbduk/CTCF_consensus_seq_D2t1r3_alpha_barcode05_sup_merge.fasta' out=/dev/null stats='/scratch/groups/astraigh/kousik/hg002_alignments/bbduk/CTCF_consensus_seq_D2t1r3_alpha_barcode05_sup_merge.txt' overwrite=t k=14 rcomp=t maskmiddle=f mink=-1 rename=t

bbduk.sh in='/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1r3_alpha/barcode06_sup/pass/barcode06_merge.fastq' -qin=33 ref='~/CTCF_consensus_seq.fasta' outm='/scratch/groups/astraigh/kousik/hg002_alignments/bbduk/CTCF_consensus_seq_D2t1r3_alpha_barcode06_sup_merge.fasta' out=/dev/null stats='/scratch/groups/astraigh/kousik/hg002_alignments/bbduk/CTCF_consensus_seq_D2t1r3_alpha_barcode06_sup_merge.txt' overwrite=t k=14 rcomp=t maskmiddle=f mink=-1 rename=t
#6/23/22

#Making backups first
mkdir /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2
mkdir /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/insitu
cd /oak/stanford/groups/astraigh/minION/all_data_backup/data

tar --use-compress-program="pigz" -cvf DiMeLo_v2/insitu/DiMeLo2_trial1_2022-06-20.tar.gz DiMeLo2_trial1 >> DiMeLo_v2/insitu/tarlog_2022-06-23.txt

ml system rclone
rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/insitu/DiMeLo2_trial1_2022-06-20.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/insitu" -v --dry-run

rclone copy /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v2/insitu/DiMeLo2_trial1_2022-06-20.tar.gz KSremote:"Straightlab Archival Data/minION_backup/DiMeLo_v2/insitu" -v

######

#6/23/22

#For fast basecalling (we use this mainly to demultiplex)

#sbatch file location: ~/sbatch_scripts/DiMeLo2_trial1_fast.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1_bar_fast.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1_bar_fast.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1_bar_fast.guppy.full.%j.err
#SBATCH --time=72:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo2_trial1/D2_trial1/20220617_1428_MN28439_FAS06729_cc3a76e7/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu' --config '/home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg' --num_callers 32 --fast5_out


#06/25/22

#first demultiplex fastq using guppy_barcoder
#submit as a sbatch

#if not sbatch
#salloc -p gpu -c 30 -G 2 --time=48:00:00
#ssh sh03-12n09 #or w.e node you get


#sbatch file location: ~/sbatch_scripts/DiMeLo2_trial1_demux.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1_bar_fast_demux.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1_bar_fast.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1_bar_fast.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g

/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_barcoder --input_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu/pass --save_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu/split --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/barcoding/configuration.cfg -t 30

sbatch ~/sbatch_scripts/DiMeLo2_trial1_demux.guppy.sbatch

---

for bar in 01 02 03 04 05 06 07 08 09 10;
do
	for file in './barcode'$bar'/*fastq';
	do
	awk -F '[ ]' 'BEGIN {OFS="\n"}(NR%4==1){print $1}' $file | sed -e 's/^.//' >> './barcode'$bar'/barcode'$bar'_readids.txt'
	done
done

'''
(base) [ksundar@sh02-ln04 login /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu/split]$ wc -l barcode*/*.txt

	158798 barcode01/barcode01_readids.txt
	345605 barcode02/barcode02_readids.txt
	331589 barcode03/barcode03_readids.txt
	442192 barcode04/barcode04_readids.txt
	370459 barcode05/barcode05_readids.txt
	266639 barcode06/barcode06_readids.txt
	295117 barcode07/barcode07_readids.txt
	226650 barcode08/barcode08_readids.txt
	238477 barcode09/barcode09_readids.txt
	97579 barcode10/barcode10_readids.txt

	2773105 total 
'''

######

##split raw fast5s by barcode readids to run through megalodon
conda activate ont2

###need to mkdir for outputdirs ./fast5_split and for each ./fast5_split/barcode0*
cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu

mkdir fast5_split

for n in 01 02 03 04 05 06 07 08 09 10;
do
mkdir 'fast5_split/barcode'$n
done


for n in 01 02 03 04 05 06 07 08 09 10;
do
fast5_subset -i /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo2_trial1/D2_trial1/20220617_1428_MN28439_FAS06729_cc3a76e7/fast5 -s '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu/fast5_split/barcode'$n -l '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu/split/barcode'$n'/barcode'$n'_readids.txt' -n 4000 -t 32 
done

'''
(ont2) [ksundar@sh02-09n13 /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu/fast5_split] (job 56361543) $ wc -l */filename_mapping.txt
	158798 barcode01/filename_mapping.txt
	345605 barcode02/filename_mapping.txt
	331589 barcode03/filename_mapping.txt
	442192 barcode04/filename_mapping.txt
	370459 barcode05/filename_mapping.txt
	266639 barcode06/filename_mapping.txt
	295117 barcode07/filename_mapping.txt
	226650 barcode08/filename_mapping.txt
	238477 barcode09/filename_mapping.txt
	97579 barcode10/filename_mapping.txt
	
	2773105 total
'''

######

###
#superaccuracy basecalling for better quality mapping

#sbatch file location: ~/sbatch_scripts/DiMeLo2_trial1_bar_supacc.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1_bar_supacc.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1_bar_supacc.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1_bar_supacc.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo2_trial1/D2_trial1/20220617_1428_MN28439_FAS06729_cc3a76e7/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1_sup/' -c dna_r9.4.1_450bps_sup.cfg -d /home/groups/astraigh/software/guppy_tar/ont-guppy/data/ --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 128 --max_queued_reads 6000 -x 'auto'


###
#for all context mod basecalling after splitting fast5s

#sbatch file location: ~/sbatch_scripts/DiMeLo2_trial1_bar_allcontext.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1_bar_allcontext.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1_bar_allcontext.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1_bar_allcontext.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo2_trial1/D2_trial1/20220617_1428_MN28439_FAS06729_cc3a76e7/fast5/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1_mod/' -c res_dna_r941_min_modbases-all-context_v001.cfg -d /home/groups/astraigh/software/rerio/basecall_models/ --bam_out --bam_methylation_threshold 0.0 --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 4096 --max_queued_reads 6000 -x 'auto'

######

#06/26/22

mkdir /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1

for n in 01 02 03 04 05 06 07 08 09 10;
do
mkdir '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/barcode'$n'_sup'
done

for n in 01 02 03 04 05 06 07 08 09 10;
do
mkdir '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/barcode'$n'_mod'
done

#sbatch scripts or superacc basecalling and all_context mod basecalling of demultiplexed libraries:

#example:
#barcode01

#Supacc basecalling

#sbatch file location: 
vim ~/sbatch_scripts/DiMeLo2_trial1_barcode01_supacc.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1_barcode01_supacc.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1_barcode01_supacc.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1_barcode01_supacc.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu/fast5_split/barcode01/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/barcode01_sup' -c dna_r9.4.1_450bps_sup.cfg -d /home/groups/astraigh/software/guppy_tar/ont-guppy/data/ --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 128 --max_queued_reads 6000 -x 'auto'

sbatch ~/sbatch_scripts/DiMeLo2_trial1_barcode01_supacc.guppy.sbatch

###
#for all context mod basecalling after splitting fast5s

#sbatch file location: 

vim ~/sbatch_scripts/DiMeLo2_trial1_barcode01_allcontext.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=DiMeLo2_trial1_barcode01_allcontext.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=DiMeLo2_trial1_barcode01_allcontext.guppy.full.%j.out
#SBATCH --error=DiMeLo2_trial1_barcode01_allcontext.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/full_fast_gpu/fast5_split/barcode01/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/barcode01_mod' -c res_dna_r941_min_modbases-all-context_v001.cfg -d /home/groups/astraigh/software/rerio/basecall_models/ --bam_out --bam_methylation_threshold 0.0 --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 4096 --max_queued_reads 6000 -x 'auto'

sbatch ~/sbatch_scripts/DiMeLo2_trial1_barcode01_allcontext.guppy.sbatch

#######

#merging modbasecalled fastq's in each barcode into 1 file

for bar in 01 02 03 04 05 06 07 08 09 10;
do
	cd '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/barcode'$bar'_mod/pass'
	for file in fastq_runid*.fastq;
	do
	cat $file >> 'barcode'$bar'_merge.fastq'
	done
done

'''
(base) [ksundar@sh02-09n13 /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1] (job 56385736) $ wc -l barcode*_mod/pass/*merge.fastq
	
	623204 barcode01_mod/pass/barcode01_merge.fastq
	1356356 barcode02_mod/pass/barcode02_merge.fastq
	1300884 barcode03_mod/pass/barcode03_merge.fastq
	1735632 barcode04_mod/pass/barcode04_merge.fastq
	1454244 barcode05_mod/pass/barcode05_merge.fastq
	1046792 barcode06_mod/pass/barcode06_merge.fastq
	1159080 barcode07_mod/pass/barcode07_merge.fastq
	890340 barcode08_mod/pass/barcode08_merge.fastq
	936468 barcode09_mod/pass/barcode09_merge.fastq
	382968 barcode10_mod/pass/barcode10_merge.fastq

	10885968 total
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

for bar in 01 02 03 04 05 06 07 08 09 10;
do
	cd '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/barcode'$bar'_sup/pass'
	for file in fastq_runid*.fastq;
	do
	cat $file >> 'barcode'$bar'_merge.fastq'
	done
done

'''
(base) [ksundar@sh02-09n13 /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1] (job 56385736) $ wc -l barcode*_sup/pass/*merge.fastq
	
	602004 barcode01_sup/pass/barcode01_merge.fastq
	1302532 barcode02_sup/pass/barcode02_merge.fastq
	1260368 barcode03_sup/pass/barcode03_merge.fastq
	1678420 barcode04_sup/pass/barcode04_merge.fastq
	1402064 barcode05_sup/pass/barcode05_merge.fastq
	1011564 barcode06_sup/pass/barcode06_merge.fastq
	1123096 barcode07_sup/pass/barcode07_merge.fastq
	860868 barcode08_sup/pass/barcode08_merge.fastq
	905812 barcode09_sup/pass/barcode09_merge.fastq
	345864 barcode10_sup/pass/barcode10_merge.fastq

	10492592 total 
'''

#winnowmap > bam; remove duplicates, unmapped etc.; sort and index

#If needed:
mkdir /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/

for bar in 01 02 03 04 05 06 07 08 09 10;
do
### KS new code: run the following
/home/groups/astraigh/software/Winnowmap/bin/winnowmap -W /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_repetitive_k15.txt -ax map-ont /oak/stanford/groups/astraigh/T2T_HG002_XY_v2.7/v2.7_HG002_w_ChrY.fasta '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/barcode'$bar'_sup/pass/barcode'$bar'_merge.fastq' | samtools view -b > '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.bam'
done

module load biology
module load samtools

for bar in 01 02 03 04 05 06 07 08 09 10;
do
#first remove unmapped (4), secondary (256), and supplemental (2048) alignments, 4+256+2048=2308
samtools view -b -F 2308 '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.bam' > /'scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.clean.bam'
#then sort and index
samtools sort '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.clean.bam' > '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.clean.sorted.bam'
samtools index '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.clean.sorted.bam'
done
######

######
#Did the following clean up immediately following winnowmap:
'''
#got to clean up the bam to avoid duplicates
module load biology
module load samtools

#first remove unmapped (4), secondary (256), and supplemental (2048) alignments, 4+256+2048=2308
for bar in 01 02 03 04 05 06 07 08 09 10;
do
samtools view -b -@ 24 -F 2308 '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.sorted.bam' > /'scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.sorted.clean.bam'
done
'''

#Making merged bam files in barcode*_mod/pass/ to sort before merging with alignment bam_out

ml biology samtools
for bar in 01 02 03 04 05 06 07 08 09 10;
do
	cd '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/barcode'$bar'_mod/pass'
	for file in *.bam;
	do
	samtools view $file >> 'barcode'$bar'_merge.txt'
	done
done

#####

##The following is OKS's code for merging the info from sup_acc + winnowmap bam and the mod_basecalled bam, modified for the current library

mkdir /scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/

for n in 01 02 03 04 05 06 07 08 09 10;
do
	mkdir '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n
	join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.9,1.10,1.11,2.12,2.13 <(samtools view -@ 30 '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/winnowmap_HG002_XY_v2.7_align.barcode'$n'_merge.clean.sorted.bam' | sort -k1) <(sort -k1 '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial1/barcode'$n'_mod/pass/barcode'$n'_merge.txt') > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.txt'
	#remove reads with empty Mm and Ml to avoid problems downstream
	grep -v 'Mm:Z:A+a;C+m;' '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.txt' > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.txt'
	#Copy header from winnowmap output bam into sam
	samtools view -H '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial1/winnowmap_HG002_XY_v2.7_align.barcode'$n'_merge.clean.sorted.bam' | cat > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.sam'
	#concatenate cleaned up merged txt file into sam file with header
	cat '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.txt' >> '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.sam'
	#convert sam to bam, filter by quality score of 10
	samtools view -@30 -b -q10 '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.sam' > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.bam'
	#sort and index bam
	samtools sort -@30 '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.bam' > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.bam'
	samtools index -@30 '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.bam'
	#generate bigwigs to look at coverage,
	/home/groups/astraigh/miniconda3/envs/deepcrap/bin/bamCoverage -p 30 -bs 100 -b '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.bam' -o '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial1/barcode'$n'/barcode'$n'_winnowmap_HG002_XY_v2.7_mod_merge.clean.filt.q10.sorted.100bp.bw'
done


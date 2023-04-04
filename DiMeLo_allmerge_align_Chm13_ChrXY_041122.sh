#Strategy for basecalling and aligning to Hybrid genome with ChrX and ChrY.

###

#From DiMeLo-seq paper:
# Basecalling for centromere enriched samples was performed twice both times using Guppy (5.0.7). The first basecalling used the “super accuracy” basecalling model (dna_r9.4.1_450bps_sup.cfg), followed by alignment to the CHM13+HG002X+hg38Y reference genome using Winnowmap (v2.03). These alignments were then filtered for only primary alignments and mapq score greater than 10 using samtools view -F 2308 -q10. A second round of basecalling was then performed again using Guppy (5.0.7) but now with the rerio all-context basecalling model (res_dna_r941_min_modbases-all-context_v001.cfg) with --bam_out and --bam_methylation_threshold 0.0. Modified basecalls were then merged by read id with winnowmap alignments to generate bam files with high confidence alignments combined with modification calls for downstream processing. For CENP-A-directed experiments four independent biological replicates were used, and for controls (IgG-directed, free-floating pA-Hia5, and untreated), two independent biological replicates were used. For all samples the first replicate was sequenced on two separate flow cells and all sequencing runs were merged for the final analysis.

###

#From OKS slack for trial12:
#For fast basecalling (we use this mainly to demultiplex)
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_1x601_730bp_trial12/trial12/20211119_1401_MN28439_FAR59876_4d4d6ba2/fast5 --save_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_1x601_730bp_trial12/full_fast_gpu --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg --num_callers 32 --fast5_out

#for all context mod basecalling after splitting fast5's
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_1x601_730bp_trial12/full_fast_gpu/fast5_split/barcode'$n --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_1x601_730bp_trial12/bar'$n'_mod' -c res_dna_r941_min_modbases-all-context_v001.cfg -d /home/groups/astraigh/software/rerio/basecall_models/ --bam_out --bam_methylation_threshold 0.0 --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 4096 --max_queued_reads 6000 -x 'auto'

#superaccuracy basecalling for better quality mapping
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_1x601_730bp_trial12/full_fast_gpu/fast5_split/barcode'$n --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_1x601_730bp_trial12/bar'$n'_sup' -c dna_r9.4.1_450bps_sup.cfg --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 128 --max_queued_reads 6000 -x 'auto'

###

##Template for using salloc and commandline 

##### guppy with basecalling and fast config for fast5 spliting 
#get resources, and ssh to compute node
salloc -p gpu -c 30 -G 2 --time=48:00:00
ssh sh03-12n09 #or w.e node you get
cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/MeAT_array_target_trial8
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path /oak/stanford/groups/astraigh/minION/all_data_backup/data/MeAT_array_target_trial8/trial8/20210422_2110_MN30093_FAP63202_981db829/fast5 --save_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/MeAT_array_target_trial8/full_fast_gpu --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg --num_callers 30 

###

##Template for sbatch instead of salloc

###RAN THIS SBATCH SCRIPT INSTEAD OF SALLOC (sbatch /home/groups/astraigh/ont_minion/scripts/sbatch/MeAT_array_trial8_guppy_full_fast_gpu.sbatch)

#!/bin/bash -l
#SBATCH --job-name=MeAT_array8.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=MeAT_array8.guppy.full.%j.out
#SBATCH --error=MeAT_array8.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path /oak/stanford/groups/astraigh/minION/all_data_backup/data/MeAT_array_target_trial8/trial8/20210422_2110_MN30093_FAP63202_981db829/fast5 --save_path /scratch/groups/astraigh/minion_seq/guppy_basecalling/MeAT_array_target_trial8/full_fast_gpu --config /home/groups/astraigh/software/guppy_tar/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg --num_callers 30    

####

#########

#04/11/2022

#Merge split fast5's of CENP-A DiMeLo together (move to same folder) from libraries below:
#Trial 1 	Barcode 18  Folderpath: "/oak/stanford/groups/astraigh/minION/DiMeLo_processed_file_backup/guppy_basecalling/DiMeLo_cen_enrich_trial1/full_fast_gpu2/fast5_split/barcode18"
#Trial 1r3 	Barcode 18 	Folderpath: "/oak/stanford/groups/astraigh/minION/DiMeLo_processed_file_backup/guppy_basecalling/DiMeLo_cen_enrich_trial1r3/full_fast_gpu2/fast5_split/barcode18"
#Trial 2 	Barcode 22 	Folderpath: "/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial2/full_fast_gpu2/fast5_split/barcode22"
#Trial 4 	Barcode 13 	Folderpath: "/oak/stanford/groups/astraigh/minION/DiMeLo_processed_file_backup/split_fast5/fast5s/DiMeLo_cen_enrich_trial4_r2/full_fast_gpu/fast5_split/barcode13.2"
#Trial 6a 	Barcode 2 	Folderpath: "/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial6a/full_fast_gpu/fast5_split/barcode02/"
#Trial 6bc 	Barcode 2 	Folderpath: "/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial6bc/full_fast_gpu/fast5_split/barcode02/"

#fast5 files are all in this folder: "/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/"

###
#Bash command for merging files:

cd /oak/stanford/groups/astraigh/minION/DiMeLo_processed_file_backup/guppy_basecalling/DiMeLo_cen_enrich_trial1/full_fast_gpu2/fast5_split/barcode18
for file in *.fast5; do cp "$file" "${file/batch/trial1_batch}"; done
mv trial1*fast5 /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/

cd /oak/stanford/groups/astraigh/minION/DiMeLo_processed_file_backup/guppy_basecalling/DiMeLo_cen_enrich_trial1r3/full_fast_gpu2/fast5_split/barcode18
for file in *.fast5; do cp "$file" "${file/batch/trial1r3_batch}"; done
mv trial1r3*fast5 /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial2/full_fast_gpu2/fast5_split/barcode22
for file in *.fast5; do cp "$file" "${file/batch/trial2_batch}"; done
mv trial2*fast5 /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/

#Permission issues for trial4, 6a, 6bc. instead moving files first then transfering:

rsync -avzh /oak/stanford/groups/astraigh/minION/DiMeLo_processed_file_backup/split_fast5/fast5s/DiMeLo_cen_enrich_trial4_r2/full_fast_gpu/fast5_split/barcode13.2/ /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial4r2 
cd /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial4r2
for file in *.fast5; do mv "$file" "${file/batch/trial4r2_batch}"; done
mv trial4r2*fast5 /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/

rsync -avzh /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial6a/full_fast_gpu/fast5_split/barcode02/ /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial6a
cd /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial6a
for file in *.fast5; do cp "$file" "${file/batch/trial6a_batch}"; done
mv trial6a*fast5 /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/

rsync -avzh /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial6bc/full_fast_gpu/fast5_split/barcode02 /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial6bc
# ^ forgot / after input folder path so have to mention barcode02 below as the folder got transferred and not just the files, minor difference
cd /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial6bc/barcode02
for file in *.fast5; do cp "$file" "${file/batch/trial6bc_batch}"; done
mv trial6bc*fast5 /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/

###filename_mapping.txt file contains readid / filename. Not sure where it is used, but making a new merged filename_mapping.txt containing readids with updated filenames

cp /oak/stanford/groups/astraigh/minION/DiMeLo_processed_file_backup/guppy_basecalling/DiMeLo_cen_enrich_trial1/full_fast_gpu2/fast5_split/barcode18/filename_mapping.txt /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial1_filename_mapping.txt
cp /oak/stanford/groups/astraigh/minION/DiMeLo_processed_file_backup/guppy_basecalling/DiMeLo_cen_enrich_trial1r3/full_fast_gpu2/fast5_split/barcode18/filename_mapping.txt /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial1r3_filename_mapping.txt
cp /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial2/full_fast_gpu2/fast5_split/barcode22/filename_mapping.txt /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial2_filename_mapping.txt
cp /oak/stanford/groups/astraigh/minION/DiMeLo_processed_file_backup/split_fast5/fast5s/DiMeLo_cen_enrich_trial4_r2/full_fast_gpu/fast5_split/barcode13.2/filename_mapping.txt /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial4r2_filename_mapping.txt
cp /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial6a/full_fast_gpu/fast5_split/barcode02/filename_mapping.txt /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial6a_filename_mapping.txt
cp /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial6bc/full_fast_gpu/fast5_split/barcode02/filename_mapping.txt /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/trial6bc_filename_mapping.txt

cd /oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/

awk '$2 = "trial1_"$2' trial1_filename_mapping.txt >> filename_mapping.txt
awk '$2 = "trial1r3_"$2' trial1r3_filename_mapping.txt >> filename_mapping.txt
awk '$2 = "trial2_"$2' trial2_filename_mapping.txt >> filename_mapping.txt
awk '$2 = "trial4r2_"$2' trial4r2_filename_mapping.txt >> filename_mapping.txt
awk '$2 = "trial6a_"$2' trial6a_filename_mapping.txt >> filename_mapping.txt
awk '$2 = "trial6bc_"$2' trial6bc_filename_mapping.txt >> filename_mapping.txt

rm -r trial*txt
#only fast5s and one filename_mapping.txt file in the folder

###
#superaccuracy basecalling for better quality mapping

#sbatch file location: ~/sbatch_scripts/CENPA_merge_supacc.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=CENPA_merge_supacc.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=CENPA_merge_supacc.guppy.full.%j.out
#SBATCH --error=CENPA_merge_supacc.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v1_CENPA_merge_sup/' -c dna_r9.4.1_450bps_sup.cfg -d /home/groups/astraigh/software/guppy_tar/ont-guppy/data/ --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 128 --max_queued_reads 6000 -x 'auto'

###
#for all context mod basecalling after splitting fast5's

#sbatch file location: ~/sbatch_scripts/CENPA_merge_allcontext.guppy.sbatch

#!/bin/bash -l
#SBATCH --job-name=CENPA_merge_allcontext.guppy.full
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ksundar@stanford.edu
#SBATCH --output=CENPA_merge_allcontext.guppy.full.%j.out
#SBATCH --error=CENPA_merge_allcontext.guppy.full.%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=gpu
#SBATCH -G 2
#SBATCH --cpus-per-task=30
##SBATCH --mem-per-cpu=15g
/home/groups/astraigh/software/guppy_tar/ont-guppy/bin/guppy_basecaller --input_path '/oak/stanford/groups/astraigh/minION/all_data_backup/data/DiMeLo_v1_CENPA_splitfast5_merge/' --save_path '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v1_CENPA_merge_mod/' -c res_dna_r941_min_modbases-all-context_v001.cfg -d /home/groups/astraigh/software/rerio/basecall_models/ --bam_out --bam_methylation_threshold 0.0 --trim_barcodes --num_callers 30 --gpu_runners_per_device 2 --chunks_per_runner 4096 --max_queued_reads 6000 -x 'auto'

######

### winnowmap to align superaccuracy fastq's to Chm13 + HG002 ChrX ChrY genome

#started doing this a while back
#2022.02.07

#chm13_T2T autosomes + MT + ChrX and ChrY v2.7 from HG002
#Downloaded latest frankengenome from https://s3-us-west-2.amazonaws.com/human-pangenomics/working/T2T/HG002XY/v2.1/v2.7.fasta locally
#Uploaded frankengenome to /oak/stanford/groups/astraigh/T2T_HG002_XY_v2.7

#preprocess kmers for T2T_HG002_XY_v2.7
cd /home/groups/astraigh/ont_minion/meryl

/home/groups/astraigh/software/Winnowmap/bin/meryl count k=15 output /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_merylDB /oak/stanford/groups/astraigh/T2T_HG002_XY_v2.7
/home/groups/astraigh/software/Winnowmap/bin/meryl print greater-than distinct=0.9998 /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_merylDB > /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_repetitive_k15.txt

######

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v1_CENPA_merge_sup/pass

for file in fastq_runid*.fastq;
do
cat $file >> CA_merge.fastq
done

wc -l CA_merge.fastq
#10574292 CA_merge.fastq

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v1_CENPA_merge_mod/pass

for file in fastq_runid*.fastq;
do
cat $file >> CA_merge_mod.fastq
done

wc -l fastq*.fastq
# 10143784
wc -l CA_merge_mod.fastq
# 10143784 CA_merge_mod.fastq

module load biology
module load samtools

### OKS code with prev version of frankengenome
#/home/groups/astraigh/software/Winnowmap/bin/winnowmap -W /home/groups/astraigh/ont_minion/meryl/chm13_HG002X_HG38Y.repetitive_k15.txt -ax map-ont /oak/stanford/groups/astraigh/charseq2.0/genomes/hsapiens/hg002_frankenstein/chm13_HG002X_HG38Y.fasta /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial2/full_fast_gpu/split/barcode23/barcode23_merge.fastq | samtools view -b > /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_cen_enrich_trial2/bar23/winnowmap_HG002align.barcode23.bam
#samtools sort winnowmap.barcode.$i.bam >winnowmap.barcode.$i.sorted.bam
#samtools index winnowmap.barcode.$i.sorted.bam
##bedtools bamtobed -i winnowmap.barcode.$i.sorted.bam > winnowmap.barcode.$i.sorted.bed
##generate bws 


### KS new code: run the following
/home/groups/astraigh/software/Winnowmap/bin/winnowmap -W /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_repetitive_k15.txt -ax map-ont /oak/stanford/groups/astraigh/T2T_HG002_XY_v2.7/v2.7_HG002_w_ChrY.fasta /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v1_CENPA_merge_sup/pass/CA_merge.fastq | samtools view -b -@ 8 > /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_v1_CA_merge/winnowmap_HG002_XY_v2.7_align.CA_merge.bam
samtools sort /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_v1_CA_merge/winnowmap_HG002_XY_v2.7_align.CA_merge.bam > /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_v1_CA_merge/winnowmap_HG002_XY_v2.7_align.CA_merge.sorted.bam
samtools index /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_v1_CA_merge/winnowmap_HG002_XY_v2.7_align.CA_merge.sorted.bam

### merge after this

join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.2,1.10,1.11,2.3,2.4 <(samtools view -@ 6 'bar'$n'_sup/winnowmap_alignHG002.v2.barcode'$n'_sup.sorted.clean.bam' | sort -k1) <(sort -k1 '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_cen_enrich_trial4_r2/bar'$n'_mod/pass/bar'$n'_mod_guppy_merge_extract.txt') > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge/DiMeLo_cen_enrich_trial4_r2/bar'$n'_sup/bar'$n'_sup_winnowmap_HG002.v2_merge.txt'

### To check coverage
module load biology
module load deeptools

bamCoverage 


#####

### kmer based classification of reads

#CENP-A enriched kmers from Logsdon et al can be found in '/home/groups/astraigh/ont_minion/unique_cen_kmer/split/chr1.chm13.v1.single.k51.mrg_sort_cenRegions.bed'

rsync -avzh ksundar@dtn.sherlock.stanford.edu:/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_v1_CA_merge/winnowmap_HG002_XY_v2.7_align.CA_merge.sorted.bam /mnt/c/Users/kousi/Desktop/Straight_lab/Straightlabcodes/Watson/alignments/dimelo

rsync -avzh ksundar@dtn.sherlock.stanford.edu:/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_v1_CA_merge/winnowmap_HG002_XY_v2.7_align.CA_merge.sorted.bam.bai /mnt/c/Users/kousi/Desktop/Straight_lab/Straightlabcodes/Watson/alignments/dimelo

###

#From Kelsey:

'''
#bw
parallel -j8 '/home/groups/astraigh/miniconda3/envs/deepcrap/bin/bamCoverage -p 30 -bs 100 -b barcode0{1}/data_splitAtME_3err_{2}.sorted.bam -o barcode0{1}/data_splitAtME_3err_{2}.sorted_100bp.bw' :::: <(seq 3 5) ::: rna dna
parallel -j8 '/home/groups/astraigh/miniconda3/envs/deepcrap/bin/bamCoverage -p 30 -bs 100 -b barcode0{1}/data_splitAtFullBridge_5err_{2}.sorted.bam -o barcode0{1}/data_splitAtFullBridge_5err_{2}.sorted_100bp.bw' :::: <(seq 3 5) ::: rna dna
#bed
parallel -j8 'bedtools bamtobed -i barcode0{1}/data_splitAtME_3err_{2}.sorted.bam > barcode0{1}/data_splitAtME_3err_{2}.sorted.bed' :::: <(seq 3 5) ::: rna dna
parallel -j8 'bedtools bamtobed -i barcode0{1}/data_splitAtFullBridge_5err_{2}.sorted.bam > barcode0{1}/data_splitAtFullBridge_5err_{2}.sorted.bed' :::: <(seq 3 5) ::: rna dna
'''

#5/24/22
#Processing 1 trial at a time (to avoid unique name issues for running annies dimelo)

#Merging fastq's

cd /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v1_CENPA_merge_sup/pass

for file in fastq_runid_1560ae9ee1d5f1206c938ef5e2e38e80b186435c*.fastq;
do
cat $file >> CA_merge_fastq_runid_1560ae9ee1d5f1206c938ef5e2e38e80b186435c.fastq
done

module load biology
module load samtools

/home/groups/astraigh/software/Winnowmap/bin/winnowmap -W /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_repetitive_k15.txt -ax map-ont /oak/stanford/groups/astraigh/T2T_HG002_XY_v2.7/v2.7_HG002_w_ChrY.fasta /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v1_CENPA_merge_sup/pass/CA_merge_fastq_runid_1560ae9ee1d5f1206c938ef5e2e38e80b186435c.fastq | samtools view -b > /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_v1_CA_merge/winnowmap_HG002_XY_v2.7_align.CA_merge_fastq_runid_1560ae9ee1d5f1206c938ef5e2e38e80b186435c.bam

'''
[M::mm_idx_gen::2.653*0.00] reading downweighted kmers                                                                                                                                                             [M::mm_idx_gen::2.760*0.02] collected downweighted kmers, no. of kmers read=66835                                                                                                                                  [M::mm_idx_gen::2.760*0.02] saved the kmers in a bloom filter: hash functions=2 and size=960936                                                                                                                    [M::mm_idx_gen::161.858*1.11] collected minimizers                                                                                                                                                                 [M::mm_idx_gen::167.623*1.15] sorted minimizers                                                                                                                                                                    [M::main::167.624*1.15] loaded/built the index for 25 target sequence(s)                                                                                                                                           [M::main::167.624*1.15] running winnowmap in SV-aware mode                                                                                                                                                         [M::main::167.624*1.15] stage1-specific parameters minP:2000, incP:2.83, maxP:16000, sample:2000, min-qlen:10000, min-qcov:0.5, min-mapq:5, mid-occ:5000                                                           [M::main::167.624*1.15] stage2-specific parameters s2_bw:2000, s2_zdropinv:25                                                                                                                                      [M::mm_idx_stat] kmer size: 15; skip: 50; is_hpc: 0; #seq: 25                                                                                                                                                      [M::mm_idx_stat::167.930*1.15] distinct minimizers: 23617830 (41.29% are singletons); average occurrences: 5.279; average spacing: 25.005                                                                          [M::worker_pipeline::8485.891*3.87] mapped 325927 sequences                                                                                                                                                        [M::worker_pipeline::16832.469*3.89] mapped 325771 sequences                                                                                                                                                       [M::worker_pipeline::25114.368*3.90] mapped 326867 sequences                                                                                                                                                       [M::worker_pipeline::31346.199*3.91] mapped 263025 sequences                                                                                                                                                       [M::main] Version: 2.03, pthreads=12, omp_threads=3                                                                                                                                                                [M::main] CMD: /home/groups/astraigh/software/Winnowmap/bin/winnowmap -W /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_repetitive_k15.txt -ax map-ont /oak/stanford/groups/astraigh/T2T_HG002_XY_v2.7/v2.7_HG002_w_ChrY.fasta /scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo_v1_CENPA_merge_sup/pass/CA_merge_fastq_runid_1560ae9ee1d5f1206c938ef5e2e38e80b186435c.fastq                           [M::main] Real time: 31348.555 sec; CPU: 122422.097 sec; Peak RSS: 46.537 GB
'''

samtools sort /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_v1_CA_merge/winnowmap_HG002_XY_v2.7_align.CA_merge_fastq_runid_1560ae9ee1d5f1206c938ef5e2e38e80b186435c.bam > /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_v1_CA_merge/winnowmap_HG002_XY_v2.7_align.CA_merge_fastq_runid_1560ae9ee1d5f1206c938ef5e2e38e80b186435c.sorted.bam

'''
[bam_sort_core] merging from 8 files and 1 in-memory blocks... 
'''

samtools index /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo_v1_CA_merge/winnowmap_HG002_XY_v2.7_align.CA_merge_fastq_runid_1560ae9ee1d5f1206c938ef5e2e38e80b186435c.sorted.bam


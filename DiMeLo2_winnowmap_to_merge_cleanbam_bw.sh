#Winnowmap_to_merge_cleanbam_bw

#Making merged bam files in barcode*_mod/pass/ to sort before merging with alignment bam_out
#Done previously

ml biology samtools
for bar in 01 02 05 06;
do
	cd '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial'$trial'_alpha/barcode'$bar'_mod/pass'
	for file in *.bam;
	do
	samtools view $file >> 'barcode'$bar'_merge.txt'
	done
done

mkdir /scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial'$trial'_alpha/
mkdir /scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/



/home/groups/astraigh/software/Winnowmap/bin/winnowmap -W /home/groups/astraigh/kousik/T2T_HG002_XY_V2.7_kmers/hg002_XY.v2.7_repetitive_k15.txt -ax map-ont /oak/stanford/groups/astraigh/T2T_HG002_XY_v2.7/v2.7_HG002_w_ChrY.fasta '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial2r1_alpha/barcode'$bar'_sup/pass/barcode'$bar'_merge.fastq' | samtools view -b > '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial2r1_alpha/winnowmap_HG002_XY_v2.7_align.barcode'$bar'_merge.bam'


#07.22.22
#This is code to be run after guppy superacc and guppy allcontext basecalling have been performed. This entire code that follows can be run as an sbatch script



module load biology
module load samtools
for trial in 1r3;
do
for bar in 01 02 05 06;
do
/home/groups/astraigh/software/Winnowmap/bin/winnowmap -W /home/groups/astraigh/kousik/T2T_HG002_verkko_v1.0_kmers/T2T_HG002_verkko_v1.0_repetitive_k15.txt -ax map-ont /oak/stanford/groups/astraigh/T2T_HG002_verkko_v1.0/T2T_HG002_verkko_v1.0.fasta '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial'$trial'_alpha/barcode'$bar'_sup/pass/barcode'$bar'_merge.fastq' | samtools view -b > '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial'$trial'_alpha/winnowmap_HG002_verkko_v1.0_align.barcode'$bar'_merge.bam'
samtools view -b -F 2308 '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial'$trial'_alpha/winnowmap_HG002_verkko_v1.0_align.barcode'$bar'_merge.bam' > /'scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial'$trial'_alpha/winnowmap_HG002_verkko_v1.0_align.barcode'$bar'_merge.clean.bam'
samtools sort '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial'$trial'_alpha/winnowmap_HG002_verkko_v1.0_align.barcode'$bar'_merge.clean.bam' > '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial'$trial'_alpha/winnowmap_HG002_verkko_v1.0_align.barcode'$bar'_merge.clean.sorted.bam'
samtools index '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial'$trial'_alpha/winnowmap_HG002_verkko_v1.0_align.barcode'$bar'_merge.clean.sorted.bam'
join -j 1 -t $'\t' -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.9,1.10,1.11,2.12,2.13 <(samtools view -@ 24 '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial'$trial'_alpha/winnowmap_HG002_verkko_v1.0_align.barcode'$n'_merge.clean.sorted.bam' | sort -k1) <(sort -k1 '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial'$trial'_alpha/barcode'$n'_mod/pass/barcode'$n'_merge.txt') > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.txt'
grep -v 'Mm:Z:A+a;C+m;' '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.txt' > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.txt'
samtools view -H '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial'$trial'_alpha/winnowmap_HG002_verkko_v1.0_align.barcode'$n'_merge.clean.sorted.bam' | cat > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.sam'
cat '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.txt' >> '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.sam'
samtools view -@30 -b -q10 '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.sam' > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.filt.q10.bam'
samtools sort -@30 '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.filt.q10.bam' > '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.filt.q10.sorted.bam'
samtools index -@30 '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.filt.q10.sorted.bam'
/home/groups/astraigh/miniconda3/envs/deepcrap/bin/bamCoverage -p 30 -bs 100 -b '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.filt.q10.sorted.bam' -o '/scratch/groups/astraigh/minion_seq/guppy_winnow_merge_KS/DiMeLo2_trial'$trial'_alpha/barcode'$n'/barcode'$n'_winnowmap_HG002_verkko_v1.0_mod_merge.clean.filt.q10.sorted.100bp.bw'
done
done



#####
for trial in 1r3;
do
for bar in 01 02 05 06;
do
/home/groups/astraigh/software/Winnowmap/bin/winnowmap -W /home/groups/astraigh/kousik/T2T_HG002_verkko_v1.0_kmers/T2T_HG002_verkko_v1.0_repetitive_k15.txt -ax map-ont /oak/stanford/groups/astraigh/T2T_HG002_verkko_v1.0/T2T_HG002_verkko_v1.0.fasta '/scratch/groups/astraigh/minion_seq/guppy_basecalling/DiMeLo2_trial'$trial'_alpha/barcode'$bar'_sup/pass/barcode'$bar'_merge.fastq' | samtools view -b > '/scratch/groups/astraigh/minion_seq/winnowmap/DiMeLo2_trial'$trial'_alpha/winnowmap_HG002_verkko_v1.0_align.barcode'$bar'_merge.bam'
done
done


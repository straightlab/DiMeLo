### READ ME ###

#DiMeLo-seq in vitro analyses
#Files in this folder were used to generate subfigures in figure 2 and extended data figures 1,2 in DiMeLo-seq manuscript (Altemose et al 2022 Nature Methods)

### Description of files:

## Figure 2 and Extended Data Figures 1,2 (in vitro DiMeLo-seq):

# DiMeLoseq_in_vitro_fastq_to_length_filtered_bam.sh

#bash script for processing fastq output of Megalodon/Guppy to make length filtered bam files. (filtering out reads > 3.6 kb ==> reads containing the full 18x601 array sequence)

# DiMeLo_1X601_730bp_CA_Array_ModProbStats_031622_Git.ipynb
# DiMeLoseq_Trial4_CA_Array_ModProbStats_2022.03.22_Git

#jupyter notebooks for (i) plotting histogram of mA/A for all reads, (ii) plotting mA/read for all positions along template, (iii) sensitivity measurements on 1x601 chromatin, (iv) clustermaps and nucleosome counting on 18x601 chromatin. Uses length filtered bam as input.
#(Trial4 - 18x601 library)

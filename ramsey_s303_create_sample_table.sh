
# list of fastqs by R1 and R2 to files from raw data folder
ls -d $PWD/*_R1*.fastq.gz > ../../../outputs/ramsey_s303/fq1.txt
ls -d $PWD/*_R2*.fastq.gz > ../../../outputs/ramsey_s303/fq2.txt

# obtain sample names from fastq paths, from both r1 and r2 pairs
cut -d'/' -f8 fq1.txt | cut -d'_' -f1,2 > ramsey_s303_sample_names.txt
cut -d'/' -f8 fq2.txt | cut -d'_' -f1,2 > ramsey_s303_sample_names_check.txt
# compare sample names from r1 and r2 to make sure correct order, no missing etc
diff ramsey_s303_sample_names.txt ramsey_s303_sample_names_check.txt

# merge sample names with paths for fastq read 1 and fastq read 2
paste ramsey_s303_sample_names.txt fq1.txt fq2.txt > ramsey_s303_sample_files.tab

# add 0s for 4th column not yet implemented for anything
awk '{print $0," 0"}' ramsey_s303_sample_files.tab > ramsey_s303_sample_table.tab

# replace tab space within "  0" with tab character made using Ctrl+V then Tab, on command line
sed -i 's/ /  /g' ramsey_s303_sample_table.tab

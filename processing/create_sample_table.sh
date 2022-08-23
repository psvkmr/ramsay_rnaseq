
# list of fastqs by R1 and R2 to files
ls -d $PWD/*R1*.fastq.gz > ../../../outputs/ramsey_s342/fq1.txt
ls -d $PWD/*R2*.fastq.gz > ../../../outputs/ramsey_s342/fq2.txt

#awk 'BEGIN { FS = "/" } ; {print $8}' fq1.txt | \
#awk 'BEGIN { FS = "_" } ; {print $1,$2,$3}' > ramsey_s342_sample_names.txt

cut -d'/' -f8 fq1.txt | cut -d'_' -f1,2 > ramsey_s342_sample_names.txt
cut -d'/' -f8 fq2.txt | cut -d'_' -f1,2 > ramsey_s342_sample_names_check.txt
diff ramsey_s342_sample_names.txt ramsey_s342_sample_names_check.txt

paste ramsey_s342_sample_names.txt fq1.txt fq2.txt > ramsey_s342_sample_files.tab

# replace tab space within "  0" with tab character made using Ctrl+V then Tab, on command line
awk '{print $0," 0"}' ramsey_s342_sample_files.tab > ramsey_s342_sample_table.tab
sed -i 's/ /  /g' ramsey_s342_sample_table.tab

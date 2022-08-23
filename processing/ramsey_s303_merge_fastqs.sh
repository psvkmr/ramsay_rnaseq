#R
# In R, read in sample table
z <- read.table('ramsey_s303_sample_table.tab')
# get unique sample IDs
uniq <- unique(z$V1)

# for each unique ID, subset the sample table, concatenate fastq file paths for each read 1
# fastq associated with that ID, store each per-sample vector of fastq file paths in vector r1
# turn into data frame
r1 <- vector()
for (n in uniq){
  v <- unique(z[z$V1 %in% n, ]$V2)
  u <- paste(v, collapse = ' ')
  r1 <- c(r1, u)
}
df1 <- data.frame(fq1 = r1)

# repeat as above with read 2 fastq files
r2 <- vector()
for (n in uniq){
  v <- unique(z[z$V1 %in% n, ]$V3)
  u <- paste(v, collapse = ' ')
  r2 <- c(r2, u)
}
df2 <- data.frame(fq2 = r2)

# write out dataframes, 1 column, 1 row per sample, each row containing all fastq files for
# read 1 or 2 for each ID
write.table(df1, 'fqs_r1.txt', row.names = F, quote = F, col.names = F)
write.table(df2, 'fqs_r2.txt', row.names = F, quote = F, col.names = F)


#shell

# create array with unique sample IDs
sample_names=ramsey_s303_sample_names.txt
unique=($(cat $sample_names | uniq))

# create array with the group of fastq files for each sample as one item
readarray fq1 < fqs_r1.txt
declare -p fq1

# create new dir to store
mkdir merge_fastqs
cd merge_fastqs

# for each sample, concatenate the multiple fastq files into 1 fastq file with all reads
# per sample for R1 reads
for i in ${!fq1[@]};
do
echo ${fq1[i]}
echo ${unique[i]}
cat ${fq1[i]} > ${unique[i]}_R1.fastq.gz
done

# repeat above process with R2 reads
cd ..

readarray fq2 < fqs_r2.txt
declare -p fq2

cd merge_fastqs

for i in ${!fq2[@]};
do
echo ${fq2[i]}
echo ${unique[i]}
cat ${fq2[i]} > ${unique[i]}_R2.fastq.gz
done

# Useful one-liners for computational biology

## Convert 

### Convert fastq to fasta

#### Option 1
```
sed -n '1~4s/^@/>/p;2~4p' file.fastq >file.fasta
```

#### Option 2
```
sed '/^@/!d;s//>/;N' file.fastq >file.fasta
```


### Convert multiline fasta to single line 

#### Option 1
```
perl -pe '/^>/ ? print "\n" : chomp' multi-line.fasta >single-line.fasta
``` 

#### Option 2
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < multi-line.fasta >single-line.fasta
``` 

#### Option 3
```
perl -pe 'chomp unless /^>/' multi-line.fasta >single-line.fasta
``` 

#### Option 4
```
awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < multi-line.fasta >single-line.fasta
```

#### Option 5
```
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' < multi-line.fasta >single-line.fasta
```

### Remove contigs based on header (after transform to single line)
```
grep -A1 -f list_to_filter contigs.fasta >rm_contigs.fasta
grep -v -f rm_contigs.fasta contigs.fasta >output.fasta
```

### Convert multi-line fasta to a single line fasta (contig files)
```
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < file.fa > out.fa
```

### Convert csv to tsv 
```
sed -E 's/("([^"]*)")?,/\2\t/g' file.csv >> file.tsv
```

### Convert Sam to fastq
```
cat file.sam | grep -v ^@ | awk '{print "@"$1"\n"$10"\n+\n"$11}' > file.fastq
```

### Copy many files from many directories
```
find . -name "*.gff" -type f -exec cp {} ./. \;
```

## Count 

### AWK for rapid fastq read count
```
awk '{s++}END{print s/4}' file.fastq
```

### Rapid fastq read count (.gz)
```
zcat file.fastq.gz | echo $((`wc -l`/4))
```

### Grep to count a fasta file

#### Option 1
```
grep ">" file.fasta -c
```

#### Option 2
```
grep ">" file.fasta | wc -l 
```

### Grep to count a fastq file

#### Option 1
```
grep "@" file.fastq -c
```

#### Option 2
```
grep "@" file.fastq | wc -l 
```

### Rename many files

#### 
```
for i in *NAME*; do mv $i ${i/NAME/NEW_NAME}; done;
```

## SED Commands

### Remove spaces in a fasta file
```
sed -r ‘s/\s+//g’
```

### Make uppercase text in a fasta file
```
sed 's/[a-z]/\U&/g'
```

### Make lowercase text in a fasta file
```
sed 's/[A-Z]/\L&/g'
``` 

### Add a tab after the first space
```
sed 's/^[ ]*\([^ ]*\) /\1\t/' file.txt >New_file.txt
```

## SAMTOOLS

### Convert Sam to Unmapped (sam)
```
samtools view -S -f 4 file.bam >file_unmapped.sam
```

### Convert Sam to Mapped (sam)
```
samtools view -S -F 4 file.sam >file_mapped.sam
```

### Convert Bam to Unmapped (bam)
```
samtools view -b -f 4 file.bam >file_unmapped.bam
```

### Convert Bam to Mapped (bam)
```
samtools view -b -F 4 file.bam >file_mapped.bam
```

### Convert Bam to fastq
```
samtools fastq file.bam >file.fastq
```

## MAPPING

## Bowtie2 

### Build Bowtie2 db
```
bowtie2-build database.fa db
```

### Map paired end/single end (.fq, .fastq, .gz)

#### Global 

##### Paired-end data (Illumina only)
```
bowtie2 -p number of threads -x db -1 R1.fastq -2 R2.fastq -S global.sam --very-sensitive
```

##### Single-end data (Illumina only)
```
bowtie2 -p number of threads -x db -q R1.fastq -S global.sam --very-sensitive
```

#### Local

##### Paired-end data (Illumina only)
```
bowtie2 -p number of threads -x db -1 1.fastq -2 2.fastq -S local.sam --very-sensitive-local
```

##### Single-end data (Illumina only)
```
bowtie2 -p number of threads -x db -q R1.fastq -S local.sam --very-sensitive-local` 
```

## Phylogenetics

### Tree annotator 
```
tree_annotator.py annotation.csv treefile output
```


### Contact 
The point-of-contact for this project is [Dr. Richard Allen White III](https://github.com/raw937).<br />
If you have any questions or feedback, please feel free to get in touch by email. 
Dr. Richard Allen White III - raw937@gmail.com.  <br />
Or [open an issue](https://github.com/raw937/Useful-one-liners-for-computational/issues).

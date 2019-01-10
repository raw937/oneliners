## Useful one-liners for computational biology

### Convert csv to tsv 
sed -E 's/("([^"]*)")?,/\2\t/g' file.csv >> file.tsv

### AWK for rapid fastq read count
awk '{s++}END{print s/4}' file.fastq

### Rapid fastq read count (.gz)
zcat file.fastq.gz | echo $((`wc -l`/4))

### Rename many files
for i in *NAME*; do mv $i ${i/NAME/NEW_NAME}; done; <br />
for i in *fastq; do mv $i ${i/R1/1}; done;

### Contact 
The point-of-contact for this project is [Dr. Richard Allen White III](https://github.com/raw937).<br />
If you have any questions or feedback, please feel free to get in touch by email. 
Dr. Richard Allen White III - raw937@gmail.com.  <br />
Or [open an issue](https://github.com/raw937/Useful-one-liners-for-computational/issues).

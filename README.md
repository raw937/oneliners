# Useful-one-liners-for-computational

##Convert csv to tsv
sed -E 's/("([^"]*)")?,/\2\t/g' file.csv >> file.tsv

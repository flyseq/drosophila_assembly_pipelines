#! /bin/bash

#set maximum number of threads to use
threads="80"

#Path to folder with all the *busco.tar.gz files
busco_path="/media/bernardkim/active-data/paper_genomes/buscos_v4_paper"

#copy all busco folders to wd and unzip
cp ${busco_path}/*.tar.gz ./
ls *.tar.gz | parallel -j${threads} tar xf {}
rm *.tar.gz

#get coordinates of complete BUSCO genes
if [ -f complete_busco_locations.csv ]; then
    rm complete_busco_locations.csv
fi
echo "Species,ID,Contig,Start,End" >> complete_busco_locations.csv
for file in $(find . -name "full_table.tsv"); do
    sp=$(echo $file | sed -E 's/\/run_diptera_odb10\/full_table.tsv//' | sed -E 's/.\///' | sed -E 's/.buscov4//')
    grep -v "^#" ${file} | awk -F'\t' -v s=${sp} '$2=="Complete" {print s","$1","$3","$4","$5}' >> complete_busco_locations.csv
done
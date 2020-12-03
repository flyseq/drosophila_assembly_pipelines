#! /bin/bash

# Script for building trees with BUSCO sequences
# This was written for doing all the analyses on a single dev box
# and makes extensive use of GNU Parallel

#set maximum number of threads to use
threads="80"

#Path to folder with all the *busco.tar.gz files
busco_path="/media/bernardkim/active-data/buscos_v4/"

#copy all busco folders to wd and unzip
#parallel -j${threads} 'pigz -dc -p 1 {} | tar xf -' ::: $(ls -lah ${busco_path}* | awk '{print $NF}' | tr '\n' ' ')

#Set the maximum number of taxa BUSCOs can be dup/missing/fragmented in
species_num=$(ls -l *.buscov4/ | grep / | wc -l)
missing_no="8"

#extract names of complete BUSCO genes
if [ -f complete_buscos.txt ]; then
    rm complete_buscos.txt
fi
for file in $(find . -name "full_table.tsv"); do
    grep -v "^#" ${file} | awk '$2=="Complete" {print $1}' >> complete_buscos.txt;
done

#filter out BUSCOs exceeding missing no threshold
sort complete_buscos.txt | uniq -c > complete_busco_ids_with_counts.txt
awk -v spnum="${species_num}" -v missing="${missing_no}" '$1 >= (spnum - missing) {print $2}' \
    complete_busco_ids_with_counts.txt > final_busco_ids.txt
echo "Number of BUSCOs used: $(cat final_busco_ids.txt | wc -l)"

#copy amino acid sequences
mkdir -p busco_aa
for dir in $(find . -type d -name "single_copy_busco_sequences"); do
  #get species names from archived busco folders
  IFS='/' read -ra a <<< "$dir"
  sppname=$(echo ${a[1]} | sed 's/.buscov4//')
  #compile amino acid info
  #to do: extract DNA sequence
  while read id; do
      file=${dir}/${id}.faa
      if [ -f $file ]; then
          #clean up fasta headers
          (echo $(echo ">${sppname}" | tr -d '.' | tr '[:lower:]' '[:upper:]' );
           cat $file | grep -v ">") > ./busco_aa/${sppname}_${id}.faa
      fi
  done < final_busco_ids.txt
done

#generate one file per BUSCO ID
while read line; do
        cat ./busco_aa/*_${line}.faa >> ./busco_aa/${line}_aa.fasta
done < final_busco_ids.txt
find . -wholename './busco_aa/*.faa' | xargs rm -f #regular rm doesn't work for many files


#align amino acids with MAFFT
mkdir -p ./busco_aa_aln
alnthreads="6" #num threads per process
parallel -j$(echo "${threads}/${alnthreads}" | bc) \
    mafft --thread ${alnthreads} --genafpair --maxiterate 10000 {} '>' \
    ./busco_aa_aln/{/.}.aln.fasta ::: $(ls -lah ./busco_aa/*.fasta | awk '{print $NF}' | tr '\n' ' ')

#Old RAxML, ignore
#mkdir raxml
#for file in ./busco_aa_aln/*.fasta; do
#    genename=(${file/aln\// })
#       genename=${genename[-1]}
#       genename=$(echo $genename | cut -f 1 -d "_")
#
#
#    raxmlHPC-PTHREADS-AVX -T 70 -f a -m PROTGAMMAAUTO -p 12345 -x 12345 -# 100 -s ${file} -n ${genename}
#done

#run raxml in parallel w/8 threads per process
#build gene trees
mkdir -p raxml
genenames=()
for file in ./busco_aa_aln/*.fasta; do
    genename=(${file/aln\// })
    genename=${genename[-1]}
    genename=$(echo $genename | cut -f 1 -d "_")
    genenames+=($genename)
done

raxthreads="5" #number of threads per process
parallel -j$(echo "${threads}/${raxthreads}"/2 | bc) raxml-ng --force --threads ${raxthreads} \
    --seed 12345 -msa ./busco_aa_aln/{}_aa.aln.fasta --msa-format FASTA --data-type AA \
    --prefix RAxML_{} --model LG '&&' mv RAxML_{}.* ./raxml/ ::: "${genenames[@]}"

#concatenate RAxML trees
cat ./raxml/*.bestTree > busco_gene_trees_ml.tree

#run astral on best trees
java -jar -D"java.library.path=lib/" ~/bin/astral.5.14.7.jar \
    --input busco_gene_trees_ml.tree \
    --output busco_species_astral.tree \
    -T 80

cat busco_species_astral.tree | sed -E 's/([A-Z\-]+[A-Z0-9]*)/\1:1/' > \
    busco_species_astral_withtiplengths.tree

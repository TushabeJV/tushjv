/home/postgrads/2488158T/kraken2-2.1.2/kraken2
/home/postgrads/2488158T/kraken2-2.1.2/kraken2-build
/home/postgrads/2488158T/kraken2-2.1.2/kraken2-inspect

wget -r -np -nH --cut-dirs=5 -R "index.html*" --accept "*.genomic.fa.gz" ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/

find kraken2_db/library/wormbase -type f -name "*.fa.gz" -exec gunzip {} +

find kraken2_db/library/wormbase -type f -name "*.fa" -exec kraken2-build --add-to-library {} --db kraken2_db \;


#****************download all to kraken3_db
*****************************************
#find kraken2_db/library/wormbase -type f -name "*.fa" -exec kraken2-build --add-to-library {} --db kraken2_db \;

#find kraken2_db/library/schistosoma -type f -name "*.fna" -exec kraken2-build --add-to-library {} --db kraken2_db \;

#find kraken2_db/library/bacteria -type f -name "*.fna" -exec kraken2-build --add-to-library {} --db kraken2_db \;

#find kraken2_db/library/plasmodium -type f -name "*.fna" -exec kraken2-build --add-to-library {} --db kraken2_db \;

#find kraken2_db/library/protozoa -type f -name "*.fna" -exec kraken2-build --add-to-library {} --db kraken2_db \;

#find kraken2_db/library/fungi -type f -name "*.fna" -exec kraken2-build --add-to-library {} --db kraken2_db \;

#find kraken2_db/library/archaea -type f -name "*.fna" -exec kraken2-build --add-to-library {} --db kraken2_db \;

#find kraken2_db/library/human -type f -name "*.fna" -exec kraken2-build --add-to-library {} --db kraken2_db \;

#find kraken2_db/library/viral -type f -name "*.fna" -exec kraken2-build --add-to-library {} --db kraken2_db \;

kraken2-build --standard --db kraken2_db --use-ftp
kraken2-build --build --db kraken2_db

kraken2_db/tar -xzvf k2_standard_8gb_20210517.tar.gz

find kraken2_db/library/bacteria/all -type f -name "*.fna.gz" -exec gunzip {} +

gunzip kraken2_db/library/bacteria/all/GCF_002947915.1_ASM294791v1_genomic.fna.gz

find kraken2_db/library/bacteria/all -type f -name "*.fna" -exec kraken2-build --add-to-library {} --db kraken2_db \;

scp -r /Users/tjv/Metagenomics/kraken2_db 2488158T@buckethead.eng.gla.ac.uk:/home/postgrads/2488158T/

/home/postgrads/2488158T/kraken2_db/library/protozoa/assembly_summary.txt": Permission denied

#!/bin/bash #PBS -N kraken2_build
# #PBS -l nodes=1:ppn=32
# #PBS -l walltime=48:00:00
# #PBS -l mem=128gb
# #PBS -o /home/postgrads/2488158T/output/kraken2_build.out
# #PBS -e /home/postgrads/2488158T/output/kraken2_build.err
# #PBS -M your.email@example.com #PBS -m abe

#***********************************************
#! /bin/bash
#$ -N split_align_merge
#$ -q postgrad.q
#$ -l mem_tokens=100.0G
#$ -o /home/postgrads/2488158T/split_align_merge.out
#$ -e /home/postgrads/2488158T/split_align_merge.out
#$ -M j.tushabe.1@research.gla.ac.uk
#$ -m beas


/home/postgrads/2488158T/kraken2-2.1.2/kraken2-build --build --db /home/postgrads/2488158T/kraken2_db --threads 32
kraken2_build.sh (END)
#******************************************

s**    /home/postgrads/2488158T/kraken2-2.1.2/kraken2-build --build --db /home/postgrads/2488158T/kraken2_db --threads 30

# moving files from computer to cluster
rsync -avz --partial --progress /Users/tjv/Metagenomics/krona_classified.fasta 2488158T@buckethead.eng.gla.ac.uk:/home/postgrads/2488158T/

# pasword for buckethead
v63ZxPhYyDf7ebgW

export LANGUAGE=en_GB.UTF-8
export LC_ALL=en_GB.UTF-8
export LC_CTYPE=en_GB.UTF-8
export LANG=en_GB.UTF-8

/home/postgrads/2488158T/kraken2-2.1.2/kraken2-inspect --db kraken2_db --threads 32 > inspectReport

# taxonomic classification
/home/postgrads/2488158T/kraken2-2.1.2/kraken2 --db kraken2_db --output kraken2_output.txt --report kraken2_report.txt --classified-out kraken2_classified.txt --unclassified-out kraken2_unclassified.txt --threads 32 --confidence 0.5 final_contigs.fa
/home/postgrads/2488158T/kraken2-2.1.2/kraken2 --db kraken2_db --output kraken2_output.txt --report kraken2_report.txt --classified-out kraken2_classified.txt --unclassified-out kraken2_unclassified.txt --threads 32 --confidence 0.5 final_contigs.fa

# moving files from cluster to my computer
rsync -avz --partial --progress 2488158T@buckethead.eng.gla.ac.uk:/home/postgrads/2488158T/krona_input.txt /Users/tjv/Metagenomics/
rsync -avz --partial --progress 2488158T@buckethead.eng.gla.ac.uk:/home/postgrads/2488158T/kraken2_report.txt /Users/tjv/Metagenomics/
rsync -avz --partial --progress 2488158T@buckethead.eng.gla.ac.uk:/home/postgrads/2488158T/id_mapping.csv /Users/tjv/Metagenomics/
rsync -avz --partial --progress 2488158T@buckethead.eng.gla.ac.uk:/home/postgrads/2488158T/kraken2_output.txt /Users/tjv/Metagenomics/
rsync -avz --partial --progress 2488158T@buckethead.eng.gla.ac.uk:/home/postgrads/2488158T/final_merged.aln /Users/tjv/Metagenomics/
rsync -avz --partial --progress 2488158T@buckethead.eng.gla.ac.uk:/home/postgrads/2488158T/split_align_merge_all.sh /Users/tjv/Metagenomics/
rsync -avz --partial --progress 2488158T@buckethead.eng.gla.ac.uk:/home/postgrads/2488158T/final_merged.ph /Users/tjv/Metagenomics/




# to run krona
~/bin/bin/ktImportTaxonomy

~/bin/bin/ktImportTaxonomy krona_input.txt -o krona_output.html

~/bin/bin/updateTaxonomy.sh

/home/postgrads/2488158T/KronaTools-2.8/updateTaxonomy.sh

rsync -avz --partial --progress 2488158T@buckethead.eng.gla.ac.uk:/home/postgrads/2488158T/krona_output.html.files /Users/tjv/Metagenomics/

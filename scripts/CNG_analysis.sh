#!bin/bash

### Define project directories
proj="50HGI"

gh_dir="$GIT_PATH/${proj}"
local_dir="$GIT_DATA/${proj}"

# Species lists based on genome quality
species_gold="${gh_dir}/auxillary/species_gold.txt"
species_ngf="${gh_dir}/auxillary/species_nongold_filarids.txt"
species_ngo="${gh_dir}/auxillary/species_nongold_other.txt"

gold_dir="${local_dir}/g_genomes"
ngf_dir="${local_dir}/ngf_genomes"
ngo_dir="${local_dir}/ngo_genomes"

### Define script output directories
# mkdir "${local_dir}/CNG/"
# mkdir "${local_dir}/CNG/g_genomes"
gold_out="${local_dir}/CNG/g_genomes"
# mkdir "${local_dir}/CNG/ngf_genomes"
ngf_out="${local_dir}/CNG/ngf_genomes"
# mkdir "${local_dir}/CNG/ngo_genomes"
ngo_out="${local_dir}/CNG/ngo_genomes"
# mkdir "${local_dir}/CNG/phylo"
phylo_out="${local_dir}/CNG/phylo"
# mkdir "${local_dir}/CNG/mcl"
# mcl_out="${local_dir}/CNG/mcl"

##Auxillary Scripts
# Extracting sequences provided list of sequence names and fasta file
seqextract_py="${gh_dir}"/scripts/auxillary/seq_extract.py
# HMMTOP parsing script (filter based on TM range and produce sequences based on TM domains)
HMMTOP_py="${gh_dir}"/scripts/auxillary/HMMTOP_extract.py
# Adding species labels to FASTA IDs
change_ID_py="${gh_dir}"/scripts/auxillary/id_change.py

## Prepare Pfam-A HMM db
# mkdir "$local_dir/auxillary"
# mkdir "$local_dir/auxillary/HMMs"
# wget -nc -O "$local_dir/auxillary/HMMs/Pfam-A.hmm.gz" ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
# gunzip -c "$local_dir/auxillary/HMMs/Pfam-A.hmm.gz" > "$local_dir/auxillary/HMMs/Pfam-A.hmm"
# hmmpress "$local_dir/auxillary/HMMs/Pfam-A.hmm"
pfam_HMM="${local_dir}"/auxillary/HMMs/Pfam-A.hmm


## CNGs
seeds="${gh_dir}"/auxillary/CNG/cng_seeds.fa

# mafft --thread 4 --auto "${seeds}" > "${gh_dir}"/auxillary/CNG/cng_seeds.aln
# trimal -gt 0.7 -cons 2 -in "${gh_dir}"/auxillary/CNG/cng_seeds.aln -out "${gh_dir}"/auxillary/CNG/cng_seeds_trim.aln
# hmmbuild "${gh_dir}"/auxillary/CNG/cng.hmm "${gh_dir}"/auxillary/CNG/cng_seeds_trim.aln
CNG_HMM="${gh_dir}"/auxillary/CNG/cng.hmm

# ### GOLD GENOMES - Prepare BLAST DBs
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 	  	#echo "${curr_dir}"
# 		gunzip -c "${f}" > "${curr_dir}"/${line}.protein.fa
# 		cd "${curr_dir}"
# 		makeblastdb -dbtype prot -in ${line}.protein.fa -out ${line}.protein.db
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - Prepare BLAST DBs
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 	  	#echo "${curr_dir}"
# 		gunzip -c "${f}" > "${curr_dir}"/${line}.protein.fa
# 		cd "${curr_dir}"
# 		makeblastdb -dbtype prot -in ${line}.protein.fa -out ${line}.protein.db
# 	done;
# done <"$species_ngf"

### GOLD GENOMES - mine for CNGs, start with reciprocal blast
# line = species name, iterate through gold genome species names
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 		# blast CNG seeds against all proteomes, using E-value cutoff
# 		cd "${curr_dir}"
# 		blastp -query "${seeds}" -db ${line}.protein.db -out "${gold_out}"/${line}.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - mine for CNGs, start with reciprocal blast
# line = species name, iterate through gold genome species names
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 		# blast CNG seeds against all proteomes, using E-value cutoff
# 		cd "${curr_dir}"
# 		blastp -query "${seeds}" -db ${line}.protein.db -out "${ngf_out}"/${line}.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species_ngf"

### GOLD GENOMES - extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cat "${gold_out}"/${line}.blastout | awk '{print $1 " " $2 " " $7}' > "${gold_out}"/${line}_blast_CNGhits.txt
# 		cat "${gold_out}"/${line}_blast_CNGhits.txt | awk '{print $2}' | awk '!seen[$0]++' > "${gold_out}"/${line}_blast_CNGhits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo ${curr_dir}
# 		gunzip -c "$f" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/${line}_blast_CNGhits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_blast_CNG.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - extract sequences of surviving hits
while IFS= read -r line; do
	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
		cat "${ngf_out}"/${line}.blastout | awk '{print $1 " " $2 " " $7}' > "${ngf_out}"/${line}_blast_CNGhits.txt
		cat "${ngf_out}"/${line}_blast_CNGhits.txt | awk '{print $2}' | awk '!seen[$0]++' > "${ngf_out}"/${line}_blast_CNGhits_ids.txt
		#Extract these sequences
		curr_dir=$(dirname "${f}")
		#echo ${curr_dir}
		gunzip -c "${f}" > "${curr_dir}"/protein.tmp.fa
		python "${seqextract_py}" "${ngf_out}"/${line}_blast_CNGhits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/${line}_blast_CNG.fa
		rm "${curr_dir}"/protein.tmp.fa
	done;
done <"$species_ngf"

### GOLD GENOMES - Reciprocal BLAST of extracted sequences against C. elegans proteome
while IFS= read -r line; do
	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
	  	curr_dir=$(dirname "${f}")
		# blast CNG seeds against all proteomes, using E-value cutoff
		cd "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/
		blastp -query "${gold_out}"/${line}_blast_CNG.fa -db caenorhabditis_elegans.protein.db -out "${gold_out}"/${line}_rec.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
	done;
done <"$species_gold"

### NON-GOLD FILARID GENOMES - Reciprocal BLAST of extracted sequences against C. elegans proteome
while IFS= read -r line; do
	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
	  	curr_dir=$(dirname "${f}")
		# blast CNG seeds against all proteomes, using E-value cutoff
		cd "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/
		blastp -query "${ngf_out}"/${line}_blast_CNG.fa -db caenorhabditis_elegans.protein.db -out "${ngf_out}"/${line}_rec.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
	done;
done <"$species_ngf"

### GOLD GENOMES - remove hits that aren't most similar to a CNG, extract sequences of surviving hits
while IFS= read -r line; do
	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
		cat "${gold_out}"/${line}_rec.blastout | awk '{print $1 " " $2 " " $7}' | sort -k1,1 -k3,3g | sort -uk1,1 | awk '/Y76B12C.1.|C23H5.7.|F38E11.12.|F14H8.6.|ZC84.2.|F36F2.5./' | sort -k3 -g > "${gold_out}"/${line}_rblast_CNGhits.txt
		cat "${gold_out}"/${line}_rblast_CNGhits.txt | awk '{print $1}' > "${gold_out}"/${line}_rblast_CNGhits_ids.txt
		#Extract these sequences
		curr_dir=$(dirname "${f}")
		#echo ${curr_dir}
		gunzip -c "${f}" > "${curr_dir}"/protein.tmp.fa
		python "${seqextract_py}" "${gold_out}"/${line}_rblast_CNGhits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_rblast_CNG.fa
		rm "${curr_dir}"/protein.tmp.fa
	done;
done <"$species_gold"

### NON-GOLD FILARID GENOMES - remove hits that aren't most similar to a CNG, extract sequences of surviving hits
while IFS= read -r line; do
	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
		cat "${ngf_out}"/${line}_rec.blastout | awk '{print $1 " " $2 " " $7}' | sort -k1,1 -k3,3g | sort -uk1,1 | awk '/Y76B12C.1.|C23H5.7.|F38E11.12.|F14H8.6.|ZC84.2.|F36F2.5./' | sort -k3 -g > "${ngf_out}"/${line}_rblast_CNGhits.txt
		cat "${ngf_out}"/${line}_rblast_CNGhits.txt | awk '{print $1}' > "${ngf_out}"/${line}_rblast_CNGhits_ids.txt
		#Extract these sequences
		curr_dir=$(dirname "${f}")
		#echo ${curr_dir}
		gunzip -c "${f}" > "${curr_dir}"/protein.tmp.fa
		python "${seqextract_py}" "${ngf_out}"/${line}_rblast_CNGhits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/${line}_rblast_CNG.fa
		rm "${curr_dir}"/protein.tmp.fa
	done;
done <"$species_ngf"

 #!bin/bash

### Define project directories
proj="Local_50HGI"

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
mkdir "${local_dir}/TRP/"
mkdir "${local_dir}/TRP/g_genomes"
gold_out="${local_dir}/TRP/g_genomes"
mkdir "${local_dir}/TRP/ngf_genomes"
ngf_out="${local_dir}/TRP/ngf_genomes"
mkdir "${local_dir}/TRP/ngo_genomes"
ngo_out="${local_dir}/TRP/ngo_genomes"
mkdir "${local_dir}/TRP/phylo"
phylo_out="${local_dir}/TRP/phylo"
mkdir "${local_dir}/TRP/mcl"
mcl_out="${local_dir}/TRP/mcl"

##Auxillary Scripts
# Extracting sequences provided list of sequence names and fasta file
seqextract_py="${gh_dir}"/scripts/auxillary/seq_extract.py
# HMMTOP parsing script (filter based on TM range and produce sequences based on TM domains)
HMMTOP_py="${gh_dir}"/scripts/auxillary/HMMTOP_extract.py
# Adding species labels to FASTA IDs
change_ID_py="${gh_dir}"/scripts/auxillary/id_change.py

### Build HMMs

## NChR s
# cat "${gh_dir}"/auxillary/pfam_HMMs/GPCR/NChR/*hmm > "${gh_dir}"/auxillary/pfam_HMMs/GPCR/NChR/NChR.hmm
# hmmpress "${gh_dir}"/auxillary/pfam_HMMs/GPCR/NChR/NChR.hmm
# NemChR_HMM="${gh_dir}"/auxillary/pfam_HMMs/GPCR/NChR/NChR.hmm

## GRAFS+
# cat "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/*hmm > "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm
# hmmpress "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm
# GRAFS_HMM="${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm

## Combined (GRAFS+ / NChR)
# cat "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm "${gh_dir}"/auxillary/pfam_HMMs/GPCR/NChR/NChR.hmm > "${gh_dir}"/auxillary/pfam_HMMs/GPCR/GRAFS_NemChR.hmm
# hmmpress "${gh_dir}"/auxillary/pfam_HMMs/GPCR/GRAFS_NemChR.hmm
# GRAFS_NemChR_HMM="${gh_dir}"/auxillary/pfam_HMMs/GPCR/GRAFS_NemChR.hmm

## Prepare Pfam-A HMM db
mkdir "$local_dir/auxillary"
mkdir "$local_dir/auxillary/HMMs"
# wget -nc -O "$local_dir/auxillary/HMMs/Pfam-A.hmm.gz" ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
# zcat "$local_dir/auxillary/HMMs/Pfam-A.hmm.gz" > "$local_dir/auxillary/HMMs/Pfam-A.hmm"
# hmmpress "$local_dir/auxillary/HMMs/Pfam-A.hmm"
pfam_HMM="${local_dir}"/auxillary/HMMs/Pfam-A.hmm

## TRPs
# seeds="${gh_dir}"/auxillary/TRP/trp_seeds.fa

# mafft --thread 4 --auto "${seeds}" > "${gh_dir}"/auxillary/TRP/trp_seeds.aln
# trimal -gt 0.7 -cons 2 -in "${gh_dir}"/auxillary/TRP/trp_seeds.aln -out "${gh_dir}"/auxillary/TRP/trp_seeds_trim.aln
# hmmbuild "${gh_dir}"/auxillary/TRP/TRP.hmm "${gh_dir}"/auxillary/TRP/trp_seeds_trim.aln
# TRP_HMM="${gh_dir}"/auxillary/TRP/TRP.hmm

### GOLD GENOMES - Prepare BLAST DBs
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 	  	#echo "${curr_dir}"
# 		zcat "${f}" > "${curr_dir}"/${line}.protein.fa
# 		cd "${curr_dir}"
# 		makeblastdb -dbtype prot -in ${line}.protein.fa -out ${line}.protein.db
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - Prepare BLAST DBs
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 	  	#echo "${curr_dir}"
# 		zcat "${f}" > "${curr_dir}"/${line}.protein.fa
# 		cd "${curr_dir}"
# 		makeblastdb -dbtype prot -in ${line}.protein.fa -out ${line}.protein.db
# 	done;
# done <"$species_ngf"

### GOLD GENOMES - mine for TRPs, start with recipricol blast
## line = species name, iterate through gold genome species names
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 		#blast TRP seeds against all proteomes, using E-value cutoff
# 		cd "${curr_dir}"
# 		blastp -query "${seeds}" -db ${line}.protein.db -out "${gold_out}"/${line}.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - mine for TRPs, start with recipricol blast
## line = species name, iterate through gold genome species names
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 		#blast TRP seeds against all proteomes, using E-value cutoff
# 		cd "${curr_dir}"
# 		blastp -query "${seeds}" -db ${line}.protein.db -out "${ngf_out}"/${line}.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species_ngf"

### GOLD GENOMES - extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cat "${gold_out}"/${line}.blastout | awk '{print $1 " " $2 " " $7}' > "${gold_out}"/${line}_blast_TRPhits.txt
# 		cat "${gold_out}"/${line}_blast_TRPhits.txt | awk '{print $2}' | awk '!seen[$0]++' > "${gold_out}"/${line}_blast_TRPhits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo ${curr_dir}
# 		zcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/${line}_blast_TRPhits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_blast_TRP.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cat "${ngf_out}"/${line}.blastout | awk '{print $1 " " $2 " " $7}' > "${ngf_out}"/${line}_blast_TRPhits.txt
# 		cat "${ngf_out}"/${line}_blast_TRPhits.txt | awk '{print $2}' | awk '!seen[$0]++' > "${ngf_out}"/${line}_blast_TRPhits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo ${curr_dir}
# 		zcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${ngf_out}"/${line}_blast_TRPhits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/${line}_blast_TRP.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_ngf"

### GOLD GENOMES - Reciprocal BLAST of extracted sequences against C. elegans proteome
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 		#blast TRP seeds against all proteomes, using E-value cutoff
# 		cd "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/
# 		blastp -query "${gold_out}"/${line}_blast_TRP.fa -db caenorhabditis_elegans.protein.db -out "${gold_out}"/${line}_rec.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - Reciprocal BLAST of extracted sequences against C. elegans proteome
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 		#blast TRP seeds against all proteomes, using E-value cutoff
# 		cd "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/
# 		blastp -query "${ngf_out}"/${line}_blast_TRP.fa -db caenorhabditis_elegans.protein.db -out "${ngf_out}"/${line}_rec.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species_ngf"

### GOLD GENOMES - remove hits that aren't most similar to a TRP, extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cat "${gold_out}"/${line}_rec.blastout | awk '{print $1 " " $2 " " $7}' | sort -k1,1 -k3,3g | sort -uk1,1 | awk '/C29E6.2.|M05B5.6.|Y73F8A.1.|ZK512.3.|T01H8.5.|F54D1.5.|C05C12.3.|ZC21.2.|R06B10.4.|K01A11.4.|R13A5.1.|B0212.5.|F28H7.10.|T09A12.3.|T10B10.7.|Y40C5A.2.|Y71A12B.4./' | sort -k3 -g > "${gold_out}"/${line}_rblast_TRPhits.txt
# 		cat "${gold_out}"/${line}_rblast_TRPhits.txt | awk '{print $1}' > "${gold_out}"/${line}_rblast_TRPhits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo ${curr_dir}
# 		zcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/${line}_rblast_TRPhits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_rblast_TRP.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - remove hits that aren't most similar to a TRP, extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cat "${ngf_out}"/${line}_rec.blastout | awk '{print $1 " " $2 " " $7}' | sort -k1,1 -k3,3g | sort -uk1,1 | awk '/C29E6.2.|M05B5.6.|Y73F8A.1.|ZK512.3.|T01H8.5.|F54D1.5.|C05C12.3.|ZC21.2.|R06B10.4.|K01A11.4.|R13A5.1.|B0212.5.|F28H7.10.|T09A12.3.|T10B10.7.|Y40C5A.2.|Y71A12B.4./' | sort -k3 -g > "${ngf_out}"/${line}_rblast_TRPhits.txt
# 		cat "${ngf_out}"/${line}_rblast_TRPhits.txt | awk '{print $1}' > "${ngf_out}"/${line}_rblast_TRPhits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo ${curr_dir}
# 		zcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${ngf_out}"/${line}_rblast_TRPhits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/${line}_rblast_TRP.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_ngf"

######
###### FIRST PHYLOGENETIC ANALYSIS
######

### GOLD GENOMES - Copy sequence files to ../phylo/ directory
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cp "${gold_out}"/${line}_rblast_TRP.fa "${phylo_out}"
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - Copy sequence files to ../phylo/ directory
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cp "${ngf_out}"/${line}_rblast_TRP.fa "${phylo_out}"
# 	done;
# done <"$species_ngf"

### Label each sequence with its species name
# for f in "${phylo_out}"/*_rblast_TRP.fa ; do
# 	python "${change_ID_py}" "$f" "$f".fa;
# done
# for f in "${phylo_out}"/*.fa.fa ; do
# 	mv "$f" "${f/.fa.fa/_label.fa}";
# done

# ### Remove Pristionchus and Panagrellus
# for f in "${phylo_out}"/panagrellus*; do mv "$f" "${f/.fa/.fa.bkp}"; done
# for f in "${phylo_out}"/pristionchus*; do mv "$f" "${f/.fa/.fa.bkp}"; done

### Cat and label files for alignment/phylo
# cat "${phylo_out}"/*_rblast_TRP_label.fa > "${phylo_out}"/All_rblast_TRP_label.fa

### Pull in manually curated outgroup (seeds)
# cat "${phylo_out}"/All_rblast_TRP_label.fa "${seeds}" > "${phylo_out}"/All_rblast_TRP_outgroup.fa


### HMMTOP
cd "${gh_dir}"/scripts/auxillary/hmmtop_2.1/
./hmmtop -if="${phylo_out}"/All_rblast_TRP_outgroup.fa -of="${phylo_out}"/All_rblast_TRP_outgroup_hmmtop_output.txt -sf=FAS
## Parse HHMTOP output to get list of seq ids with >= 5 TM domains <= 10 TM domains
~/install/anaconda3/bin/python "${HMMTOP_py}" "${phylo_out}"/All_NCf_outgroup_hmmtop_output.txt "${phylo_out}"/All_NCf_outgroup.fa "${phylo_out}"/All_NCf_outgroup_TMfiltered.fa
## Remove sequences without any predicted TMs
cat "${phylo_out}"/All_rblast_TRP_outgroup_hmmtop_output.txt | awk '{print $3 " " $5}' | awk '$2!=0' | awk '{print $1}' > "${phylo_out}"/All_rblast_TRP_outgroup.key
~/install/anaconda3/bin/python "${seqextract_py}" "${phylo_out}"/All_rblast_TRP_outgroup.key "${phylo_out}"/All_rblast_TRP_outgroup.fa "${phylo_out}"/All_rblast_TRPf_outgroup.fa


### Align files
einsi --thread 4 "${phylo_out}"/All_rblast_TRPf_outgroup.fa > "${phylo_out}"/All_rblast_TRPf_outgroup.aln

### Trim alignments
# trimal_cmd="${gh_dir}"/scripts/auxillary/trimal/source/./trimal
# "${trimal_cmd}" -in "${phylo_out}"/All_rblast_TRPf_outgroup.aln -out "${phylo_out}"/All_rblast_TRPf_outgroup_trim.aln -gt 0.75 -cons 2

# on server
# ssh zamanian@brc6.secure.biotech.wisc.edu 'nohup ~/install/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -f a -x 12345 -p 12345 -# 100 -m PROTCATAUTO -s ~/data/TRP/All_rblast_TRPf_outgroup_trim.aln -n TRP &'


######
###### MCL CLUSTERING ANALYSIS
######

# cp "${phylo_out}"/All_TRPf_outgroup.fa "${mcl_out}"

# make_db="${gh_dir}"/scripts/auxillary/makeblastdb
# blast="${gh_dir}/"scripts/auxillary/blastp
# load="${gh_dir}/"scripts/auxillary/mcxload
# mcl="${gh_dir}/"scripts/auxillary/mcl
# mcxdump="${gh_dir}/"scripts/auxillary/mcxdump
# prepCSV_py="${gh_dir}/"scripts/auxillary/format_csv.py

## blast all vs all
# cd "${mcl_out}"
# "${make_db}" -dbtype prot -in All_TRPf_outgroup.fa -out All_TRPf_outgroup  
# "${blast}" -db All_TRPf_outgroup -query All_TRPf_outgroup.fa -out blastall.out -evalue 0.01 -outfmt 6
## prepare for MCL
# cut -f 1,2,11 blastall.out > blastall.abc
# "${load}" -abc blastall.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o blastall.mci -write-tab blastall.tab
## MCL analysis
# "${mcl}" blastall.mci -I 1.2
# "${mcl}" blastall.mci -I 1.4
# "${mcl}" blastall.mci -I 2
# "${mcl}" blastall.mci -I 4
# "${mcl}" blastall.mci -I 6

# "${mcxdump}" -icl out.blastall.mci.I12 -tabr blastall.tab -o dump.blastall.mci.I12
# "${mcxdump}" -icl out.blastall.mci.I14 -tabr blastall.tab -o dump.blastall.mci.I14
# "${mcxdump}" -icl out.blastall.mci.I20 -tabr blastall.tab -o dump.blastall.mci.I20
# "${mcxdump}" -icl out.blastall.mci.I40 -tabr blastall.tab -o dump.blastall.mci.I40
# "${mcxdump}" -icl out.blastall.mci.I60 -tabr blastall.tab -o dump.blastall.mci.I60

# python "${prepCSV_py}" "${mcl_out}"/dump.blastall.mci.I12 "${mcl_out}"/clusters.csv












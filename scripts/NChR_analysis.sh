#!bin/bash

### Define project directories
proj="50HGI"

gh_dir="${GIT_PATH}"/"${proj}"
local_dir="${GIT_DATA}"/"${proj}"

# Species lists based on genome quality
species_gold="${gh_dir}"/auxillary/species_gold.txt
species_ngf="${gh_dir}"/auxillary/species_nongold_filarids.txt
species_ngo="${gh_dir}"/auxillary/species_nongold_other.txt
chemoR_families="${gh_dir}"/auxillary/chemoR_families.txt

gold_dir="${local_dir}"/g_genomes
ngf_dir="${local_dir}"/ngf_genomes
ngo_dir="${local_dir}"/ngo_genomes

### Define script output directories
# mkdir "${local_dir}"/NChR/
# mkdir "${local_dir}"/NChR/g_genomes
# mkdir "${local_dir}"/NChR/ngf_genomes
# mkdir "${local_dir}"/NChR/ngo_genomes
# mkdir "${local_dir}"/NChR/phylo
# mkdir "${local_dir}"/NChR/mcl

gold_out="${local_dir}"/NChR/g_genomes
ngo_out="${local_dir}"/NChR/ngo_genomes
ngf_out="${local_dir}"/NChR/ngf_genomes
phylo_out="${local_dir}"/NChR/phylo
mcl_out="${local_dir}"/NChR/mcl

##Auxillary Scripts
# Extracting sequences provided list of sequence names and fasta file
seqextract_py="${gh_dir}"/scripts/auxillary/seq_extract.py
# HMMTOP parsing script (filter based on TM range and produce sequences based on TM domains)
HMMTOP_py="${gh_dir}"/scripts/auxillary/HMMTOP_extract.py
HMMTOP_strict_py="${gh_dir}"/scripts/auxillary/HMMTOP_extract_strict.py
# Adding species labels to FASTA IDs
change_ID_py="${gh_dir}"/scripts/auxillary/id_change.py
trimal_cmd="${gh_dir}"/scripts/auxillary/trimal/source/./trimal
mafft_cmd="${gh_dir}"/scripts/auxillary/mafft
muscle_cmd="${gh_dir}"/scripts/auxillary/muscle
makeblastdb_cmd="${gh_dir}"/scripts/auxillary/makeblastdb
blastp_cmd="${gh_dir}"/scripts/auxillary/blastp

### Build HMMs

## NChRs
# cat "${gh_dir}"/auxillary/pfam_HMMs/GPCR/NChR/*hmm > "${gh_dir}"/auxillary/pfam_HMMs/GPCR/NChR/NChR.hmm
# hmmpress "${gh_dir}"/auxillary/pfam_HMMs/GPCR/NChR/NChR.hmm
NemChR_HMM="${gh_dir}"/auxillary/pfam_HMMs/GPCR/NChR/NChR.hmm

## GRAFS+
# cat "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/*hmm > "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm
# hmmpress "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm
GRAFS_HMM="${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm

## Combined (GRAFS+ / NChR)
# cat "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm "${gh_dir}"/auxillary/pfam_HMMs/GPCR/NChR/NChR.hmm > "${gh_dir}"/auxillary/pfam_HMMs/GPCR/GRAFS_NemChR.hmm
# hmmpress "${gh_dir}"/auxillary/pfam_HMMs/GPCR/GRAFS_NemChR.hmm
GRAFS_NemChR_HMM="${gh_dir}"/auxillary/pfam_HMMs/GPCR/GRAFS_NemChR.hmm

# Prepare Pfam-A HMM db
# mkdir "$local_dir/auxillary
# mkdir "$local_dir/auxillary/HMMs
# wget -nc -O "$local_dir/auxillary/HMMs/Pfam-A.hmm.gz" ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.g
# gzcat "$local_dir/auxillary/HMMs/Pfam-A.hmm.gz" > "$local_dir/auxillary/HMMs/Pfam-A.hmm
# hmmpress "$local_dir/auxillary/HMMs/Pfam-A.hmm
pfam_HMM="${local_dir}"/auxillary/HMMs/Pfam-A.hmm


################################################################################################################################
###########																											 ###########
###########																											 ###########
###########										ChemoR Identification												 ###########
###########																											 ###########
###########																											 ###########
################################################################################################################################


### GOLD GENOMES - mine for nematode chemo Rs
## line = species name, iterate through gold genome species names
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 	  	#echo "${curr_dir}"
# 		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa

# 		#HMMSEARCH all proteomes against db of All GPCR hmms
# 		hmmsearch --tblout "${gold_out}"/"${line}"_hits.out --noali "${GRAFS_NemChR_HMM}" "${curr_dir}"/protein.tmp.fa 

# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

# ## NON-GOLD FILARID GENOMES - mine for nematode chemo Rs
# line = species name, iterate through gold genome species names
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 	  	#echo "${curr_dir}"
# 		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa

# 		#HMMSEARCH all proteomes against db of All GPCR hmms
# 		hmmsearch --tblout "${ngf_out}"/"${line}"_hits.out --noali "${GRAFS_NemChR_HMM}" "${curr_dir}"/protein.tmp.fa 

# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_ngf"


### GOLD GENOMES - Parse hmm outputs to filter out those where first hit is not NChR hmm, extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cat "${gold_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_1|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g > "${gold_out}"/"${line}"_NChits.txt
# 		cat "${gold_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_1|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g | awk '{print $1}' > "${gold_out}"/"${line}"_NChits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo "${curr_dir}"
# 		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/"${line}"_NChits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/"${line}"_NC.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

# ## NON-GOLD FILARID GENOMES - Parse hmm outputs to filter out those where first hit is not NChR hmm, extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cat "${ngf_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_1|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g > "${ngf_out}"/"${line}"_NChits.txt
# 		cat "${ngf_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_1|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g | awk '{print $1}' > "${ngf_out}"/"${line}"_NChits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo "${curr_dir}"
# 		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${ngf_out}"/"${line}"_NChits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/"${line}"_NC.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_ngf"

### GOLD GENOMES - Reciprocal HMMSEARCH of extracted sequences against pfam-a
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		#HMMSEARCH all putative chemosensory genes against db of PFAM hmms
# 		hmmsearch --tblout "${gold_out}"/"${line}"_rHMM.out --noali "${pfam_HMM}" "${gold_out}"/"${line}"_NC.fa
# 	done;
# done <"$species_gold"

# ## NON-GOLD FILARID GENOMES - Reciprocal HMMSEARCH of extracted sequences against pfam-a
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		#HMMSEARCH all putative chemosensory genes against db of PFAM hmms
# 		hmmsearch --tblout "${ngf_out}"/"${line}"_rHMM.out --noali "${pfam_HMM}" "${ngf_out}"/"${line}"_NC.fa
# 	done;
# done <"$species_ngf"
		

### GOLD GENOMES - Parse hmm outputs to remove sequences where first hit is not NChR, get list of surviving unique seq ids, extract sequences 
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cat "${gold_out}"/"${line}"_rHMM.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g > "${gold_out}"/"${line}"_NChitsf.txt
# 		cat "${gold_out}"/"${line}"_rHMM.out | awk '{print $1 " " $3  " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g  | awk '{print $1}' > "${gold_out}"/"${line}"_NChitsf_ids.txt
# 		curr_dir=$(dirname "${f}")
#  		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/"${line}"_NChitsf_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/"${line}"_NCf.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 		grep -v -f "${gold_out}"/"${line}"_NChitsf_ids.txt "${gold_out}"/"${line}"_NChits_ids.txt > "${gold_out}"/"${line}"_filtered1_ids.tx
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - Parse hmm outputs to remove sequences where first hit is not NChR, get list of surviving unique seq ids, extract sequences 
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cat "${ngf_out}"/"${line}"_rHMM.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g > "${ngf_out}"/"${line}"_NChitsf.txt
# 		cat "${ngf_out}"/"${line}"_rHMM.out | awk '{print $1 " " $3  " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g  | awk '{print $1}' > "${ngf_out}"/"${line}"_NChitsf_ids.txt
# 		curr_dir=$(dirname "${f}")
#  		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${ngf_out}"/"${line}"_NChitsf_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/"${line}"_NCf.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 		grep -v -f "${ngf_out}"/"${line}"_NChitsf_ids.txt "${ngf_out}"/"${line}"_NChits_ids.txt > "${ngf_out}"/"${line}"_filtered1_ids.tx
# 	done;
# done <"$species_ngf"

### GOLD GENOMES - Parse hmm outputs for each species to get list of unique seq ids for R-A, R-P, G, F, A/S (can also use  if ($4 <= 1e-50) in awk)
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		#Rhodopsin
# 		cat "${gold_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g > "${gold_out}"/"${line}"_Rhits.txt
# 		cat "${gold_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g  | awk '{if ($4 <= 1e-0) print $1}' > "${gold_out}"/"${line}"_Rhits_ids.txt
# 		curr_dir=$(dirname "${f}")
#  		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/"${line}"_Rhits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/"${line}"_R.fa
# 		#Glutamate
# 		cat "${gold_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_3/' | sort -k4 -g > "${gold_out}"/"${line}"_Ghits.txt
# 		cat "${gold_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_3/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${gold_out}"/"${line}"_Ghits_ids.txt
# 		python "${seqextract_py}" "${gold_out}"/"${line}"_Ghits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/"${line}"_G.fa
# 		#Adhesion/Secretin
# 		cat "${gold_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_2/' | sort -k4 -g > "${gold_out}"/"${line}"_AShits.txt
# 		cat "${gold_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_2/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${gold_out}"/"${line}"_AShits_ids.txt
# 		python "${seqextract_py}" "${gold_out}"/"${line}"_AShits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/"${line}"_AS.fa
# 		#Frizzled
# 		cat "${gold_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/Frizzled/' | sort -k4 -g > "${gold_out}"/"${line}"_Fhits.txt
# 		cat "${gold_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/Frizzled/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${gold_out}"/"${line}"_Fhits_ids.txt
# 		python "${seqextract_py}" "${gold_out}"/"${line}"_Fhits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/"${line}"_F.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - Parse hmm outputs for each species to get list of unique seq ids for R-A, R-P, G, F, A/S (can also use  if ($4 <= 1e-50) in awk)
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		#Rhodopsin
# 		cat "${ngf_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g > "${ngf_out}"/"${line}"_Rhits.txt
# 		cat "${ngf_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g  | awk '{if ($4 <= 1e-0) print $1}' > "${ngf_out}"/"${line}"_Rhits_ids.txt
# 		curr_dir=$(dirname "${f}")
#  		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${ngf_out}"/"${line}"_Rhits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/"${line}"_R.fa
# 		#Glutamate
# 		cat "${ngf_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_3/' | sort -k4 -g > "${ngf_out}"/"${line}"_Ghits.txt
# 		cat "${ngf_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_3/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${ngf_out}"/"${line}"_Ghits_ids.txt
# 		python "${seqextract_py}" "${ngf_out}"/"${line}"_Ghits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/"${line}"_G.fa
# 		#Adhesion/Secretin
# 		cat "${ngf_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_2/' | sort -k4 -g > "${ngf_out}"/"${line}"_AShits.txt
# 		cat "${ngf_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_2/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${ngf_out}"/"${line}"_AShits_ids.txt
# 		python "${seqextract_py}" "${ngf_out}"/"${line}"_AShits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/"${line}"_AS.fa
# 		#Frizzled
# 		cat "${ngf_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/Frizzled/' | sort -k4 -g > "${ngf_out}"/"${line}"_Fhits.txt
# 		cat "${ngf_out}"/"${line}"_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/Frizzled/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${ngf_out}"/"${line}"_Fhits_ids.txt
# 		python "${seqextract_py}" "${ngf_out}"/"${line}"_Fhits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/"${line}"_F.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_ngf"

### GOLD GENOMES - Reciprocal blastp of extracted sequences against C. elegans
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		#blast filtered NChRs against C. elegans proteome, using E-value cutoff
# 		cd "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/
# 		"${blastp_cmd}" -query "${gold_out}"/"${line}"_NCf.fa -db caenorhabditis_elegans.protein.db -out "${gold_out}"/"${line}"_rec.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - Reciprocal blastp of extracted sequences against C. elegans
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		#blast filtered NChRs against C. elegans proteome, using E-value cutoff
# 		cd "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/
# 		"${blastp_cmd}" -query "${ngf_out}"/"${line}"_NCf.fa -db caenorhabditis_elegans.protein.db -out "${ngf_out}"/"${line}"_rec.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species_ngf"


### GOLD GENOMES - remove hits that aren't most similar to a C. elegans ChemoR, extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cat "${gold_out}"/"${line}"_rec.blastout | awk '{print $1 " " $2 " " $7}' | sort -k1,1 -k3,3g | sort -uk1,1 | grep -wF -f "${gold_out}"/caenorhabditis_elegans_NChitsf_ids.txt | sort -k3 -g > "${gold_out}"/"${line}"_rblast_ChemoRhits.txt
# 		cat "${gold_out}"/"${line}"_rblast_ChemoRhits.txt | awk '{print $1}' > "${gold_out}"/"${line}"_rblast_ChemoRhits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo "${curr_dir}"
# 		zcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/"${line}"_rblast_ChemoRhits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/"${line}"_rblast_ChemoR.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 		grep -v -f "${gold_out}"/"${line}"_rblast_ChemoRhits_ids.txt "${gold_out}"/"${line}"_NChitsf_ids.txt > "${gold_out}"/"${line}"_filtered2_ids.txt
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - remove hits that aren't most similar to a C. elegans ChemoR, extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cat "${ngf_out}"/"${line}"_rec.blastout | awk '{print $1 " " $2 " " $7}' | sort -k1,1 -k3,3g | sort -uk1,1 | grep -wF -f "${gold_out}"/caenorhabditis_elegans_NChitsf_ids.txt | sort -k3 -g > "${ngf_out}"/"${line}"_rblast_ChemoRhits.txt
# 		cat "${ngf_out}"/"${line}"_rblast_ChemoRhits.txt | awk '{print $1}' > "${ngf_out}"/"${line}"_rblast_ChemoRhits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo "${curr_dir}"
# 		zcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${ngf_out}"/"${line}"_rblast_ChemoRhits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/"${line}"_rblast_ChemoR.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 		grep -v -f "${ngf_out}"/"${line}"_rblast_ChemoRhits_ids.txt "${ngf_out}"/"${line}"_NChitsf_ids.txt > "${ngf_out}"/"${line}"_filtered2_ids.txt
# 	done;
# done <"$species_ngf"

### Double check filtered2 for filarids (filtered2_description.xlsx). Retain 8 of the 105 that were filtered; manually add to the ChemoRhits IDs and FASTA file

################################################################################################################################
###########																											 ###########
###########																											 ###########
###########												Phylogenetics												 ###########
###########																											 ###########
###########																											 ###########
################################################################################################################################

### GOLD GENOMES - Copy sequence files to ../phylo/NemChR directory
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cp "${gold_out}"/"${line}"_rblast_ChemoR.fa "${phylo_out}"
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - Copy sequence files to ../phylo/NemChR directory
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cp "${ngf_out}"/"${line}"_rblast_ChemoR.fa "${phylo_out}"
# 	done;
# done <"$species_ngf"

### Label each sequence with its species name
# for f in "${phylo_out}"/*_rblast_ChemoR.fa ; do
# 	python "${change_ID_py}" "$f" "$f".fa;
# done

# for f in "${phylo_out}"/*.fa.fa ; do
# 	mv "$f" "${f/.fa.fa/_label.fa}";
# done

### Choose one or more representatives from each non-filarid clade and concatenate (can add outgroups if wanted)
# cat "${phylo_out}"/caenorhabditis_elegans_rblast_ChemoR_label.fa "${phylo_out}"/necator_americanus_rblast_ChemoR_label.fa "${phylo_out}"/haemonchus_contortus_rblast_ChemoR_label.fa "${phylo_out}"/strongyloides_ratti_rblast_ChemoR_label.fa "${phylo_out}"/trichinella_spiralis_rblast_ChemoR_label.fa "${phylo_out}"/toxocara_canis_rblast_ChemoR_label.fa > "${phylo_out}"/DS_non-filarid.fa
# cat "${phylo_out}"/brugia_pahangi_rblast_ChemoR_label.fa "${phylo_out}"/wuchereria_bancrofti_rblast_ChemoR_label.fa "${phylo_out}"/onchocerca_ochengi_rblast_ChemoR_label.fa "${phylo_out}"/brugia_timori_rblast_ChemoR_label.fa "${phylo_out}"/dirofilaria_immitis_rblast_ChemoR_label.fa  "${phylo_out}"/brugia_malayi_rblast_ChemoR_label.fa "${phylo_out}"/loa_loa_rblast_ChemoR_label.fa "${phylo_out}"/onchocerca_volvulus_rblast_ChemoR_label.fa  > "${phylo_out}"/DS_filarid.fa

# Join files
# cat "${phylo_out}"/DS_filarid.fa "${phylo_out}"/DS_non-filarid.fa > "${phylo_out}"/DS_NC.fa

## HMMTOP
# cd "${gh_dir}"/scripts/auxillary/hmmtop_2.1/
# ./hmmtop -if="${phylo_out}"/DS_non-filarid.fa -of="${phylo_out}"/DS_non-filarid_hmmtop_output.txt -sf=FAS
# ./hmmtop -if="${phylo_out}"/DS_filarid.fa -of="${phylo_out}"/DS_filarid_hmmtop_output.txt -sf=FAS

### Parse HHMTOP output to get FASTA file of non-filarid sequences with exactly 7 TMs; extract entire sequence,; don't filter filarids
# python "${HMMTOP_strict_py}" "${phylo_out}"/DS_non-filarid_hmmtop_output.txt "${phylo_out}"/DS_non-filarid.fa "${phylo_out}"/DS_non-filarid_TMfiltered.fa

### prepare C. elegans ChemoR fasta file from ___________
# awk '/^>/ {printf("%s%s\n",(N>0?"\n":"),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < "${local_dir}"/auxillary/12915_2008_195_MOESM33_ESM.fast" > "${phylo_out}"/BMC.celeg.fasta

### Create a separate fasta file for each family
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/sra/' > "${phylo_out}"/sra.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srab/' >"${phylo_out}"/srab.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srh/' > "${phylo_out}"/srh.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/str/' > "${phylo_out}"/str.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/sri/' > "${phylo_out}"/sri.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srd/' > "${phylo_out}"/srd.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srj/' > "${phylo_out}"/srj.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srm/' > "${phylo_out}"/srm.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srn/' > "${phylo_out}"/srn.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/sre/' > "${phylo_out}"/sre.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srb/' > "${phylo_out}"/srb.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srx/' > "${phylo_out}"/srx.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srt/' > "${phylo_out}"/srt.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srg/' > "${phylo_out}"/srg.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/sru/' > "${phylo_out}"/sru.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srv/' > "${phylo_out}"/srv.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srxa/' >"${phylo_out}"/srxa.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srw/' > "${phylo_out}"/srw.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srz/' > "${phylo_out}"/srz.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srbc/' >"${phylo_out}"/srbc.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srsx/' >"${phylo_out}"/srsx.celeg.fa
# cat "${phylo_out}"/BMC.celeg.fasta | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/srr/' > "${phylo_out}"/srr.celeg.fa

# for f in "${phylo_out}"/s*.celeg.fa; do
# 	cat "${f}" | tr '\t' '\n' > $f.fasta;
# 	done
# for f in "${phylo_out}"/*.fa.fasta ; do
# 	mv "$f" "${f/.fa.fasta/.fa}";
# 	done

### align family fasta files
# for f in "${phylo_out}"/s*.celeg.fa; do "${mafft_cmd}" --auto "$f" > "$f.aln"; done
# for f in "${phylo_out}"/*.fa.aln ; do
# 	mv "$f" "${f/.fa.aln/.aln}";
# 	done

### align family profiles to each other, sequentially
# muscle -profile -in1 "${phylo_out}"/sra.celeg.aln -in2 "${phylo_out}"/srab.celeg.aln -out "${phylo_out}"/2.aln
# muscle -profile -in1 "${phylo_out}"/2.aln -in2 "${phylo_out}"/srh.celeg.aln -out "${phylo_out}"/3.aln
# muscle -profile -in1 "${phylo_out}"/3.aln -in2 "${phylo_out}"/str.celeg.aln -out "${phylo_out}"/4.aln
# muscle -profile -in1 "${phylo_out}"/4.aln -in2 "${phylo_out}"/sri.celeg.aln -out "${phylo_out}"/5.aln
# muscle -profile -in1 "${phylo_out}"/5.aln -in2 "${phylo_out}"/srd.celeg.aln -out "${phylo_out}"/6.aln
# muscle -profile -in1 "${phylo_out}"/6.aln -in2 "${phylo_out}"/srj.celeg.aln -out "${phylo_out}"/7.aln
# muscle -profile -in1 "${phylo_out}"/7.aln -in2 "${phylo_out}"/sre.celeg.aln -out "${phylo_out}"/8.aln
# muscle -profile -in1 "${phylo_out}"/8.aln -in2 "${phylo_out}"/srb.celeg.aln -out "${phylo_out}"/9.aln
# muscle -profile -in1 "${phylo_out}"/9.aln -in2 "${phylo_out}"/srx.celeg.aln -out "${phylo_out}"/10.aln
# muscle -profile -in1 "${phylo_out}"/10.aln -in2 "${phylo_out}"/srt.celeg.aln -out "${phylo_out}"/11.aln
# muscle -profile -in1 "${phylo_out}"/11.aln -in2 "${phylo_out}"/srg.celeg.aln -out "${phylo_out}"/12.aln
# muscle -profile -in1 "${phylo_out}"/12.aln -in2 "${phylo_out}"/sru.celeg.aln -out "${phylo_out}"/13.aln
# muscle -profile -in1 "${phylo_out}"/13.aln -in2 "${phylo_out}"/srxa.celeg.aln -out "${phylo_out}"/14.aln
# muscle -profile -in1 "${phylo_out}"/14.aln -in2 "${phylo_out}"/srw.celeg.aln -out "${phylo_out}"/15.aln
# muscle -profile -in1 "${phylo_out}"/15.aln -in2 "${phylo_out}"/srz.celeg.aln -out "${phylo_out}"/16.aln
# muscle -profile -in1 "${phylo_out}"/16.aln -in2 "${phylo_out}"/srbc.celeg.aln -out "${phylo_out}"/17.aln
# muscle -profile -in1 "${phylo_out}"/17.aln -in2 "${phylo_out}"/srsx.celeg.aln -out "${phylo_out}"/18.aln
# muscle -profile -in1 "${phylo_out}"/18.aln -in2 "${phylo_out}"/srv.celeg.aln -out "${phylo_out}"/19.aln

#manually remove all celeg sequences from DS_non-filarid_TMfiltered2.fa

### add non-C.elegans representatives to the alignment
# "${mafft_cmd}" --reorder --thread 8 --addfull "${phylo_out}"/DS_non-filarid_TMfiltered2.fa --keeplength "${phylo_out}"/19.aln > "${phylo_out}"/non-filarid.aln
### add filarid ChemoRs to the alignment
# "${mafft_cmd}" --reorder --thread 8 --addfull "${phylo_out}"/DS_filarid.fa --keeplength "${phylo_out}"/non-filarid.aln > "${phylo_out}"/final.aln
### trim and filter
# "${trimal_cmd}" -gt 0.7 -in "${phylo_out}"/final.aln -out "${phylo_out}"/final.trim.aln
# "${trimal_cmd}" -resoverlap 0.70 -seqoverlap 70 -in "${phylo_out}"/final.trim.aln -out "${phylo_out}"/final.trim.filter.aln
### Change to single-line FASTA
# awk '/^>/ {printf("%s%s\n",(N>0?"\n":"),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < "${phylo_out}"/final.trim.aln > "${phylo_out}"/final.trim-single.aln
# awk '/^>/ {printf("%s%s\n",(N>0?"\n":"),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < "${phylo_out}"/final.trim.filter.aln | sed '/^$/d' > "${phylo_out}"/final.trim.filter-single.aln
### Get IDs and compare lists
# cat "${phylo_out}"/final.trim-single.aln | tr '\t' '\n' | awk 'NR%2==1' > "${phylo_out}"/trimmed_ids.txt
# cat "${phylo_out}"/final.trim.filter-single.aln | tr '\t' '\n' | awk 'NR%2==1' > "${phylo_out}"/filtered_ids.txt
# grep -v -f "${phylo_out}"/filtered_ids.txt "${phylo_out}"/trimmed_ids.txt > "${phylo_out}"/missing_ids.txt
# cat "${phylo_out}"/final.trim.filter.aln | sed 's/ce|/celeg-/' > "${phylo_out}"/final.trim.filter2.align
# mv "${phylo_out}"/final.trim.filter2.align "${phylo_out}"/final.trim.filter.aln

### ML tree on server
# /home/BIOTECH/zamanian/install/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 4 -f a -x 12345 -p 12345 -# 100 -m PROTGAMMAVT -s final.trim.filter.aln -n final_ML


################################################################################################################################
###########																											 ###########
###########																											 ###########
###########												Heatmap														 ###########
###########																											 ###########
###########																											 ###########
################################################################################################################################


# python "${seqextract_py}" "${phylo_out}"/ML/family_clades_ids.txt "${local_dir}"/p_genomes/phylo.fa "${phylo_out}"/ML/family_clades.fa
## manually edit to remove protein names (eg. srw-25) and species abbreviations
# auxillary/makeblastdb  -dbtype prot -in "${phylo_out}"/ML/family_clades.fa -out "${phylo_out}"/ML/family_clades.db

### Label each sequence with its species name
# for f in "${gold_out}"/*_rblast_ChemoR.fa ; do
# 	python "${change_ID_py}" "$f" "$f".fa;
# done

# for f in "${gold_out}"/*.fa.fa ; do
# 	mv "$f" "${f/.fa.fa/_label.fa}";
# done

### Label each sequence with its species name
# for f in "${ngf_out}"/*_rblast_ChemoR.fa ; do
# 	python "${change_ID_py}" "$f" "$f".fa;
# done

# for f in "${ngf_out}"/*.fa.fa ; do
# 	mv "$f" "${f/.fa.fa/_label.fa}";
# done


### GOLD GENOMES - BLAST species not included in original tree against families
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cd "${phylo_out}"/ML/
# 		"${blastp_cmd}" -query "${gold_out}"/"${line}"_rblast_ChemoR_label.fa -db family_clades.db -out "${gold_out}"/"${line}"_ChemoR_family.blastout -num_threads 4 -evalue 0.01 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID - BLAST species not included in original tree against families
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cd "${phylo_out}"/ML/
# 		"${blastp_cmd}" -query "${ngf_out}"/"${line}"_rblast_ChemoR_label.fa -db family_clades.db -out "${ngf_out}"/"${line}"_ChemoR_family.blastout -num_threads 4 -evalue 0.01 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species_ngf"

## concatenate blast results and remove self to self HSPs
# cat "${gold_out}"/*_ChemoR_family.blastout "${ngf_out}"/*_ChemoR_family.blastout | sed 's/-/\t/' | awk '$2!=$3' | sed 's/lloa_/lloa/g' > "${phylo_out}"/ChemoR_parse.blastout

################################################################################################################################
###########																											 ###########
###########																											 ###########
###########											Family Alignments												 ###########
###########																											 ###########
###########																											 ###########
################################################################################################################################

## Label each sequence with its species name
# for f in "${local_dir}"/a_genomes/*.fa; do
# 	python "${change_ID_py}" "$f" "$f".fa
# done

# for f in "${local_dir}"/a_genomes/*.fa.fa; do
# 	mv "$f" "${f/.fa.fa/_label.fa}"
# done

# cat "${local_dir}"/a_genomes/*_label.fa > "${local_dir}"/a_genomes/all.fa

# while IFS= read -r line; do
# 	python "${seqextract_py}" "${phylo_out}"/ML/families/"${line}"_ids.txt "${local_dir}"/a_genomes/all.fa "${phylo_out}"/ML/families/"${line}".fa
# 	"${mafft_cmd}" --thread 4 --reorder --auto "${phylo_out}"/ML/families/"${line}".fa > "${phylo_out}"/ML/families/"${line}".aln
# 	"${trimal_cmd}" -gt 0.7 -in "${phylo_out}"/ML/families/"${line}".aln -out "${phylo_out}"/ML/families/"${line}".trim.aln
# 	"${trimal_cmd}" -resoverlap 0.70 -seqoverlap 70 -in "${phylo_out}"/ML/families/"${line}".trim.aln -out "${phylo_out}"/ML/families/"${line}".trim.filter.aln
# done <"$chemoR_families"

# align family profiles to each other to make superfamily alignments
# NOTE: OS X overwrites files with same name but different capitalization
# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/srh.trim.aln -in2 "${phylo_out}"/ML/families/str.trim.aln -out "${phylo_out}"/ML/families/Str1.aln
# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/Str1.aln -in2 "${phylo_out}"/ML/families/sri.trim.aln -out "${phylo_out}"/ML/families/Str2.aln
# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/Str2.aln -in2 "${phylo_out}"/ML/families/srd.trim.aln -out "${phylo_out}"/ML/families/Str3.aln
# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/Str3.aln -in2 "${phylo_out}"/ML/families/srj.trim.aln -out "${phylo_out}"/ML/families/Str_sf.aln

# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/sre.trim.aln -in2 "${phylo_out}"/ML/families/sra.trim.aln -out "${phylo_out}"/ML/families/Sra1.aln
# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/Sra1.aln -in2 "${phylo_out}"/ML/families/srab.trim.aln -out "${phylo_out}"/ML/families/Sra2.aln
# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/Sra2.aln -in2 "${phylo_out}"/ML/families/srb.trim.aln -out "${phylo_out}"/ML/families/Sra_sf.aln

# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/srx.trim.aln -in2 "${phylo_out}"/ML/families/srt.trim.aln -out "${phylo_out}"/ML/families/Srg1.aln
# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/Srg1.aln -in2 "${phylo_out}"/ML/families/srg.trim.aln -out "${phylo_out}"/ML/families/Srg2.aln
# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/Srg2.aln -in2 "${phylo_out}"/ML/families/sru.trim.aln -out "${phylo_out}"/ML/families/Srg3.aln
# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/Srg3.aln -in2 "${phylo_out}"/ML/families/srv.trim.aln -out "${phylo_out}"/ML/families/Srg4.aln
# "${muscle_cmd}" -profile -in1 "${phylo_out}"/ML/families/Srg4.aln -in2 "${phylo_out}"/ML/families/srxa.trim.aln -out "${phylo_out}"/ML/families/Srg_sf.aln

# rm "${phylo_out}"/ML/families/S*1.aln
# rm "${phylo_out}"/ML/families/S*2.aln
# rm "${phylo_out}"/ML/families/S*3.aln
# rm "${phylo_out}"/ML/families/S*4.aln

while IFS= read -r line; do
	cd "${phylo_out}"/ML/families/
	/home/BIOTECH/zamanian/install/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 4 -f a -x 12345 -p 12345 -# 100 -m PROTGAMMAVT -s "${line}.trim.filter.aln" -n "${line}".tre;
done <"$chemoR_families"

# while IFS= read -r line; do
# 	awk '/^>/ {printf("%s%s\n",(N>0?"\n":"),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < "${phylo_out}"/ML/families/"${line}".trim.aln > "${phylo_out}"/ML/families/"${line}".trim-single.aln
# 	awk '/^>/ {printf("%s%s\n",(N>0?"\n":"),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < "${phylo_out}"/ML/families/"${line}".trim.filter.aln > "${phylo_out}"/ML/families/"${line}".trim.filter-single.aln
# 	cat "${phylo_out}"/ML/families/"${line}".trim-single.aln | tr '\t' '\n' | awk 'NR%2==1' > "${phylo_out}"/ML/families/"${line}".trimmed_ids.txt
# 	cat "${phylo_out}"/ML/families/"${line}".trim.filter-single.aln | tr '\t' '\n' | awk 'NR%2==1' > "${phylo_out}"/ML/families/"${line}".filtered_ids.txt
# 	grep -v -f "${phylo_out}"/ML/families/"${line}".filtered_ids.txt "${phylo_out}"/ML/families/"${line}".trimmed_ids.txt > "${phylo_out}"/ML/families/"${line}".missing_ids.txt
# done <"$chemoR_families"

# rm "${phylo_out}"/ML/families/*single*

# while IFS= read -r line; do
# 	cd "${phylo_out}"/ML/families
# 	cat "${phylo_out}"/ML/families/"${line}".trim.filter.aln | sed 's/-//' | "${makeblastdb_cmd}" -dbtype prot -in - -out "${line}".db -title "${line}"
# done <"$chemoR_families"

# while IFS= read -r line; do
# 	cat "${phylo_out}"/ML/families/"${line}".missing_ids.txt | sed 's/>//' | python "${seqextract_py}" /dev/stdin "${local_dir}"/a_genomes/all.fa "${phylo_out}"/ML/families/"${line}".missing.fa
# 	cd "${phylo_out}"/ML/families/
# 	"${blastp_cmd}" -query "${line}".missing.fa -db "${line}".db -out "${line}".missing.blastout -num_threads 4 -evalue 0.01 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# done <"$chemoR_families"
































#!bin/bash

### Define project directories
proj="50HGI"

gh_dir=${GIT_PATH}/${proj}
local_dir=${GIT_DATA}/${proj}

# Species lists based on genome quality
species_gold="${gh_dir}/auxillary/species_gold.txt"
species_ngf="${gh_dir}/auxillary/species_nongold_filarids.txt"
species_ngo="${gh_dir}/auxillary/species_nongold_other.txt"

gold_dir="${local_dir}/g_genomes"
ngf_dir="${local_dir}/ngf_genomes"
ngo_dir="${local_dir}/ngo_genomes"

### Define script output directories
# mkdir "${local_dir}/NChR/"
# mkdir "${local_dir}/NChR/g_genomes"
# mkdir "${local_dir}/NChR/ngf_genomes"
# mkdir "${local_dir}/NChR/ngo_genomes"
# mkdir "${local_dir}/NChR/phylo"
# mkdir "${local_dir}/NChR/mcl"

gold_out="${local_dir}/NChR/g_genomes"
ngo_out="${local_dir}/NChR/ngo_genomes"
ngf_out="${local_dir}/NChR/ngf_genomes"
phylo_out="${local_dir}/NChR/phylo"
mcl_out="${local_dir}/NChR/mcl"

##Auxillary Scripts
# Extracting sequences provided list of sequence names and fasta file
seqextract_py="${gh_dir}"/scripts/auxillary/seq_extract.py
# HMMTOP parsing script (filter based on TM range and produce sequences based on TM domains)
HMMTOP_py="${gh_dir}"/scripts/auxillary/HMMTOP_extract.py
HMMTOP_strict_py="${gh_dir}"/scripts/auxillary/HMMTOP_extract_strict.py
# Adding species labels to FASTA IDs
change_ID_py="${gh_dir}"/scripts/auxillary/id_change.py

### Build HMMs

## NChR s
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
# mkdir "$local_dir/auxillary"
# mkdir "$local_dir/auxillary/HMMs"
# wget -nc -O "$local_dir/auxillary/HMMs/Pfam-A.hmm.gz" ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
# gzcat "$local_dir/auxillary/HMMs/Pfam-A.hmm.gz" > "$local_dir/auxillary/HMMs/Pfam-A.hmm"
# hmmpress "$local_dir/auxillary/HMMs/Pfam-A.hmm"
pfam_HMM="$local_dir/auxillary/HMMs/Pfam-A.hmm"


### GOLD GENOMES - mine for nematode chemo Rs
## line = species name, iterate through gold genome species names
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 	  	#echo "${curr_dir}"
# 		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa

# 		#HMMSEARCH all proteomes against db of All GPCR hmms
# 		hmmsearch --tblout "${gold_out}"/${line}_hits.out --noali "${GRAFS_NemChR_HMM}" "${curr_dir}"/protein.tmp.fa 

# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

# ## NON-GOLD FILARID GENOMES - mine for nematode chemo Rs
# line = species name, iterate through gold genome species names
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 	  	#echo "${curr_dir}"
# 		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa

# 		#HMMSEARCH all proteomes against db of All GPCR hmms
# 		hmmsearch --tblout "${ngf_out}"/${line}_hits.out --noali "${GRAFS_NemChR_HMM}" "${curr_dir}"/protein.tmp.fa 

# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_ngf"


### GOLD GENOMES - Parse hmm outputs to filter out those where first hit is not NChR hmm, extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cat "${gold_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_1|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g > "${gold_out}"/${line}_NChits.txt
# 		cat "${gold_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_1|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g | awk '{print $1}' > "${gold_out}"/${line}_NChits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo ${curr_dir}
# 		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/${line}_NChits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_NC.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

# ## NON-GOLD FILARID GENOMES - Parse hmm outputs to filter out those where first hit is not NChR hmm, extract sequences of surviving hits
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cat "${ngf_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_1|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g > "${ngf_out}"/${line}_NChits.txt
# 		cat "${ngf_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_1|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g | awk '{print $1}' > "${ngf_out}"/${line}_NChits_ids.txt
# 		#Extract these sequences
# 		curr_dir=$(dirname "${f}")
# 		#echo ${curr_dir}
# 		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${ngf_out}"/${line}_NChits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/${line}_NC.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_ngf"

### GOLD GENOMES - Reciprocal HMMSEARCH of extracted sequences against pfam-a
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		#HMMSEARCH all putative chemosensory genes against db of PFAM hmms
# 		hmmsearch --tblout "${gold_out}"/${line}_rHMM.out --noali "${pfam_HMM}" "${gold_out}"/${line}_NC.fa
# 	done;
# done <"$species_gold"

# ## NON-GOLD FILARID GENOMES - Reciprocal HMMSEARCH of extracted sequences against pfam-a
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 		#HMMSEARCH all putative chemosensory genes against db of PFAM hmms
# 		hmmsearch --tblout "${ngf_out}"/${line}_rHMM.out --noali "${pfam_HMM}" "${ngf_out}"/${line}_NC.fa
# 	done;
# done <"$species_ngf"
		

### GOLD GENOMES - Parse hmm outputs to remove sequences where first hit is not NChR, get list of surviving unique seq ids, extract sequences 
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cat "${gold_out}"/${line}_rHMM.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g > "${gold_out}"/${line}_NChitsf.txt
# 		cat "${gold_out}"/${line}_rHMM.out | awk '{print $1 " " $3  " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g  | awk '{print $1}' > "${gold_out}"/${line}_NChitsf_ids.txt
# 		curr_dir=$(dirname "${f}")
#  		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/${line}_NChitsf_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_NCf.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

# ### NON-GOLD FILARID GENOMES - Parse hmm outputs to remove sequences where first hit is not NChR, get list of surviving unique seq ids, extract sequences 
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cat "${ngf_out}"/${line}_rHMM.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g > "${ngf_out}"/${line}_NChitsf.txt
# 		cat "${ngf_out}"/${line}_rHMM.out | awk '{print $1 " " $3  " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g  | awk '{print $1}' > "${ngf_out}"/${line}_NChitsf_ids.txt
# 		curr_dir=$(dirname "${f}")
#  		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${ngf_out}"/${line}_NChitsf_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/${line}_NCf.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_ngf"

##### Remove C. elegans pseudogenes (not working)
### rm "${gold_out}"/caenorhabditis_elegans_NCf.fa
##
#### extract pseudogenes from C. elegans GTF
### gtf_parser="${gh_dir}"/scripts/auxillary/gtf_parse.py
### gzcat "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS9.canonical_geneset.gtf.gz > "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS9.canonical_geneset.gtf
### python "${gtf_parser}" "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS9.canonical_geneset.gtf > "${local_dir}"/auxillary/GTF_pseudogenes.txt
### awk '{print $6}' "${local_dir}"/auxillary/GTF_pseudogenes.txt > "${local_dir}"/auxillary/GTF_pseudogenes_Tid.txt
### rm "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS9.canonical_geneset.gtf
#### compare list of pseudogenes to list of NCf hits; remove pseudogenes from list of NCf
### grep -Fxv -f "${local_dir}"/auxillary/GTF_pseudogenes_Tid.txt "${gold_out}"/caenorhabditis_elegans_NChitsf_ids.txt > "${gold_out}"/caenorhabditis_elegans_NChitsff_ids.txt

### GOLD GENOMES - Parse hmm outputs for each species to get list of unique seq ids for R-A, R-P, G, F, A/S (can also use  if ($4 <= 1e-50) in awk)
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		#Rhodopsin
# 		cat "${gold_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g > "${gold_out}"/${line}_Rhits.txt
# 		cat "${gold_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g  | awk '{if ($4 <= 1e-0) print $1}' > "${gold_out}"/${line}_Rhits_ids.txt
# 		curr_dir=$(dirname "${f}")
#  		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/${line}_Rhits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_R.fa
# 		#Glutamate
# 		cat "${gold_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_3/' | sort -k4 -g > "${gold_out}"/${line}_Ghits.txt
# 		cat "${gold_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_3/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${gold_out}"/${line}_Ghits_ids.txt
# 		python "${seqextract_py}" "${gold_out}"/${line}_Ghits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_G.fa
# 		#Adhesion/Secretin
# 		cat "${gold_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_2/' | sort -k4 -g > "${gold_out}"/${line}_AShits.txt
# 		cat "${gold_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_2/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${gold_out}"/${line}_AShits_ids.txt
# 		python "${seqextract_py}" "${gold_out}"/${line}_AShits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_AS.fa
# 		#Frizzled
# 		cat "${gold_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/Frizzled/' | sort -k4 -g > "${gold_out}"/${line}_Fhits.txt
# 		cat "${gold_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/Frizzled/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${gold_out}"/${line}_Fhits_ids.txt
# 		python "${seqextract_py}" "${gold_out}"/${line}_Fhits_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_F.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - Parse hmm outputs for each species to get list of unique seq ids for R-A, R-P, G, F, A/S (can also use  if ($4 <= 1e-50) in awk)
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 		#Rhodopsin
# 		cat "${ngf_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g > "${ngf_out}"/${line}_Rhits.txt
# 		cat "${ngf_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g  | awk '{if ($4 <= 1e-0) print $1}' > "${ngf_out}"/${line}_Rhits_ids.txt
# 		curr_dir=$(dirname "${f}")
#  		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${ngf_out}"/${line}_Rhits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/${line}_R.fa
# 		#Glutamate
# 		cat "${ngf_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_3/' | sort -k4 -g > "${ngf_out}"/${line}_Ghits.txt
# 		cat "${ngf_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_3/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${ngf_out}"/${line}_Ghits_ids.txt
# 		python "${seqextract_py}" "${ngf_out}"/${line}_Ghits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/${line}_G.fa
# 		#Adhesion/Secretin
# 		cat "${ngf_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_2/' | sort -k4 -g > "${ngf_out}"/${line}_AShits.txt
# 		cat "${ngf_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_2/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${ngf_out}"/${line}_AShits_ids.txt
# 		python "${seqextract_py}" "${ngf_out}"/${line}_AShits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/${line}_AS.fa
# 		#Frizzled
# 		cat "${ngf_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/Frizzled/' | sort -k4 -g > "${ngf_out}"/${line}_Fhits.txt
# 		cat "${ngf_out}"/${line}_hits.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/Frizzled/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${ngf_out}"/${line}_Fhits_ids.txt
# 		python "${seqextract_py}" "${ngf_out}"/${line}_Fhits_ids.txt "${curr_dir}"/protein.tmp.fa "${ngf_out}"/${line}_F.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_ngf"

make_db="${gh_dir}"/scripts/auxillary/makeblastdb
# "${make_db}" -dbtype prot -in "${gold_dir}/caenorhabditis_elegans/PRJNA13758"/caenorhabditis_elegans.protein.fa -out "${gold_dir}/caenorhabditis_elegans/PRJNA13758"/caenorhabditis_elegans.protein

### GOLD GENOMES - Reciprocal blastp of extracted sequences against C. elegans
while IFS= read -r line; do
	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
		#blast filtered NChRs against C. elegans proteome, using E-value cutoff
		cd "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/
		blastp -query "${gold_out}"/${line}_NCf.fa -db caenorhabditis_elegans.protein.db -out "${gold_out}"/${line}_rec.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
	done;
done <"$species_gold"

## NON-GOLD FILARID GENOMES - Reciprocal blastp of extracted sequences against C. elegans
while IFS= read -r line; do
	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
		#blast filtered NChRs against C. elegans proteome, using E-value cutoff
		cd "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/
		blastp -query "${ngf_out}"/${line}_NCf.fa -db caenorhabditis_elegans.protein.db -out "${ngf_out}"/${line}_rec.blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
	done;
done <"$species_ngf"


######
###### PHYLOGENETIC ANALYSIS
######

### GOLD GENOMES - Copy sequence files to ../phylo/NemChR directory
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cp "${gold_out}"/${line}_*.fa "${phylo_out}"
# 	done;
# done <"$species_gold"

### NON-GOLD FILARID GENOMES - Copy sequence files to ../phylo/NemChR directory
# while IFS= read -r line; do
# 	for f in "${ngf_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cp "${ngf_out}"/${line}_*.fa "${phylo_out}"
# 	done;
# done <"$species_ngf"

### Remove the unfiltered NC files
# rm "${phylo_out}"/*_NC.fa

### Label each sequence with its species name
# for f in "${phylo_out}"/*.fa ; do
# 	python "${change_ID_py}" "$f" "$f".fa;
# done
# for f in "${phylo_out}"/*.fa.fa ; do
# 	mv "$f" "${f/.fa.fa/_label.fa}";
# done

### Remove Pristionchus and Panagrellus
# for f in "${phylo_out}"/panagrellus*; do mv "$f" "${f/.fa/.fa.bkp}"; done
# for f in "${phylo_out}"/pristionchus*; do mv "$f" "${f/.fa/.fa.bkp}"; done

### Cat and label files for alignment/phylo
# cat "${phylo_out}"/*_NCf_label.fa > "${phylo_out}"/All_NCf.fa
# cat "${phylo_out}"/*_R_label.fa > "${phylo_out}"/All_R.fa
# cat "${phylo_out}"/*_G_label.fa > "${phylo_out}"/All_G.fa
# cat "${phylo_out}"/*_F_label.fa > "${phylo_out}"/All_F.fa
# cat "${phylo_out}"/*_AS_label.fa > "${phylo_out}"/All_AS.fa

### Label each sequence with its predicted family
# sed 's/>/>NC:/' "${phylo_out}"/All_NCf.fa > "${phylo_out}"/All_NCf_lab.fa
# sed 's/>/>R:/' "${phylo_out}"/All_R.fa > "${phylo_out}"/All_Rf_lab.fa
# sed 's/>/>G:/' "${phylo_out}"/All_G.fa > "${phylo_out}"/All_Gf_lab.fa
# sed 's/>/>F:/' "${phylo_out}"/All_F.fa > "${phylo_out}"/All_Ff_lab.fa
# sed 's/>/>AS:/' "${phylo_out}"/All_AS.fa > "${phylo_out}"/All_ASf_lab.fa
# cat "${phylo_out}"/All*_lab.fa > "${phylo_out}"/All.fa

### Pull in manually curated outgroup
# outgroup_fa="${gh_dir}/auxillary/NChR/outgroup.fa"
# cat "${phylo_out}"/All_NCf.fa "${outgroup_fa}" > "${phylo_out}"/All_NCf_outgroup.fa


### HMMTOP
# cd "${gh_dir}"/scripts/auxillary/hmmtop_2.1/
# ./hmmtop -if="${phylo_out}"/All_NCf_outgroup.fa -of="${phylo_out}"/All_NCf_outgroup_hmmtop_output.txt -sf=FAS

## Parse HHMTOP output to get list of seq ids with >= 5 TM domains <= 10 TM domains
# python "${HMMTOP_py}" "${phylo_out}"/All_NCf_outgroup_hmmtop_output.txt "${phylo_out}"/All_NCf_outgroup.fa "${phylo_out}"/All_NCf_outgroup_TMfiltered.fa

### Align files
#mafft --op 2 --ep 1 --thread 2 --maxiterate 1 "${phylo_out}"/All_NCf_outgroup_TMfiltered.fa > "${phylo_out}"/All_NCf_outgroup_TMfiltered_align.fa
#mafft --op 2 --ep 1 --thread 2 --maxiterate 1 "${phylo_out}"/All_NCf_outgroup.fa > "${phylo_out}"/All_NCf_outgroup_align.fa

### Trim alignments
trimal_cmd="${gh_dir}"/scripts/auxillary/trimal/source/./trimal
# "${trimal_cmd}" -in "${phylo_out}"/All_NCf_outgroup_TMfiltered.aln -out "${phylo_out}"/All_NCf_outgroup_TMfiltered_trim.aln -gt 0.8 -cons 2
# "${trimal_cmd}" -in "${phylo_out}"/All_NCf_outgroup.aln -out "${phylo_out}"/All_NCf_outgroup_trim.aln -gt 0.8 -cons 2


######
###### MCL CLUSTERING ANALYSIS
######

### Re-Add Pristionchus and Panagrellus
# for f in "${phylo_out}"/panagrellus*; do mv "$f" "${f/.fa.bkp/.fa}"; done
# for f in "${phylo_out}"/pristionchus*; do mv "$f" "${f/.fa.bkp/.fa}"; done

# cat "${phylo_out}"/*_NCf_label.fa "${outgroup_fa}" > "${mcl_out}"/All_NCf_outgroup.fa

# make_db="${gh_dir}"/scripts/auxillary/makeblastdb
# blast="${gh_dir}/"scripts/auxillary/blastp
# load="${gh_dir}/"scripts/auxillary/mcxload
# mcl="${gh_dir}/"scripts/auxillary/mcl
# mcxdump="${gh_dir}/"scripts/auxillary/mcxdump
# prepCSV_py="${gh_dir}/"scripts/auxillary/format_csv.py

## blast all vs all
# cd "${mcl_out}"
# "${make_db}" -dbtype prot -in All_NCf_outgroup.fa -out All_NCf_outgroup  
# "${blast}" -db All_NCf_outgroup -query All_NCf_outgroup.fa -out blastall.out -evalue 0.01 -outfmt 6
## prepare for MCL
# cut -f 1,2,11 blastall.out > blastall.abc
# "${load}" -abc blastall.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o blastall.mci -write-tab blastall.tab
## MCL analysis
# "${mcl}" blastall.mci -I 1.4
# "${mcl}" blastall.mci -I 2
# "${mcl}" blastall.mci -I 4
# "${mcl}" blastall.mci -I 6
# "${mcl}" blastall.mci -I 1.2

# "${mcxdump}" -icl out.blastall.mci.I14 -tabr blastall.tab -o dump.blastall.mci.I14
# "${mcxdump}" -icl out.blastall.mci.I20 -tabr blastall.tab -o dump.blastall.mci.I20
# "${mcxdump}" -icl out.blastall.mci.I40 -tabr blastall.tab -o dump.blastall.mci.I40
# "${mcxdump}" -icl out.blastall.mci.I60 -tabr blastall.tab -o dump.blastall.mci.I60
# "${mcxdump}" -icl out.blastall.mci.I12 -tabr blastall.tab -o dump.blastall.mci.I12

# python "${prepCSV_py}" "${mcl_out}"/dump.blastall.mci.I12 "${mcl_out}"/clusters.csv



###
### DOWN-SAMPLED PHYLOGENETIC TREE
###


### Choose one or more representatives from each non-filarid clade and concatenate (can add outgroups if wanted)
# cat "${phylo_out}"/caenorhabditis_elegans_NCf_label.fa "${phylo_out}"/necator_americanus_NCf_label.fa "${phylo_out}"/haemonchus_contortus_NCf_label.fa "${phylo_out}"/strongyloides_ratti_NCf_label.fa "${phylo_out}"/trichinella_spiralis_NCf_label.fa "${phylo_out}"/toxocara_canis_NCf_label.fa > "${phylo_out}"/DS_non-filarid.fa
# cat "${phylo_out}"/brugia_pahangi_NCf_label.fa "${phylo_out}"/wuchereria_bancrofti_NCf_label.fa "${phylo_out}"/onchocerca_ochengi_NCf_label.fa "${phylo_out}"/brugia_timori_NCf_label.fa "${phylo_out}"/dirofilaria_immitis_NCf_label.fa  "${phylo_out}"/brugia_malayi_NCf_label.fa "${phylo_out}"/loa_loa_NCf_label.fa "${phylo_out}"/onchocerca_volvulus_NCf_label.fa  > "${phylo_out}"/DS_filarid.fa
# sed 's/>/>out-/' "${outgroup_fa}" > "${phylo_out}"/outgroup.fa
# cat "${phylo_out}"/DS_non-filarid.fa "${phylo_out}"/outgroup.fa > "${phylo_out}"/DS_non-filarid_outgroup.fa

## HMMTOP
# cd "${gh_dir}"/scripts/auxillary/hmmtop_2.1/
# ./hmmtop -if="${phylo_out}"/DS_non-filarid.fa -of="${phylo_out}"/DS_non-filarid_hmmtop_output.txt -sf=FAS
# ./hmmtop -if="${phylo_out}"/DS_filarid.fa -of="${phylo_out}"/DS_filarid_hmmtop_output.txt -sf=FAS

### Parse HHMTOP output to get list of seq ids with >= 3 TM domains <= 10 TM domains for filarids, and 7 TM domains for non-filarids
# extract only TM domains
# python "${HMMTOP_strict_py}" "${phylo_out}"/DS_non-filarid_hmmtop_output.txt "${phylo_out}"/DS_non-filarid.fa "${phylo_out}"/DS_non-filarid_TMfiltered.fa
# python "${HMMTOP_py}" "${phylo_out}"/DS_filarid_hmmtop_output.txt "${phylo_out}"/DS_filarid.fa "${phylo_out}"/DS_filarid_TMfiltered.fa

# Join files
# cat "${phylo_out}"/DS_filarid_TMfiltered.fa "${phylo_out}"/DS_non-filarid_TMfiltered.fa > "${phylo_out}"/DS_NC2.fa

### Align files
# einsi --thread 4 "${phylo_out}"/DS_NC2.fa > "${phylo_out}"/DS_NC2.aln
### Email when complete
# mailx -s "Alignment complete!" njwheeler@wisc.edu <<< "The alignment of DS_NC2 has successfully completed."

### Trim alignment
## fix loa loa IDs
# sed 's/lloa_/lloa/' "${phylo_out}"/DS_NC2.aln > "${phylo_out}"/DS_NC3.aln
# mv "${phylo_out}"/DS_NC3.aln "${phylo_out}"/DS_NC2.aln
## Trim
# trimal_cmd="${gh_dir}"/scripts/auxillary/trimal/source/./trimal
# "${trimal_cmd}" -in "${phylo_out}"/DS_NC2.aln -out "${phylo_out}"/DS_NC2_trim.aln -gt 0.7
## Split in to DS_NC2_filarid_trim.aln and DS_NC2_non-filarid_trim.aln, to remove short poorly aligned sequences from non-filarids
# cat "${phylo_out}"/DS_NC2_trim.aln | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '/bmala|bpaha|wbanc|ooche|btimo|dimmi|lloa|ovolv/' > "${phylo_out}"/DS_NC2_filarid_trim.aln
# cat "${phylo_out}"/DS_NC2_trim.aln | awk '/>.*$/ { printf("%s\t", $0); next } 1' | awk '!/bmala|bpaha|wbanc|ooche|btimo|dimmi|lloa|ovolv/' > "${phylo_out}"/DS_NC2_non-filarid_trim.aln
# cat "${phylo_out}"/DS_NC2_filarid_trim.aln | tr '\t' '\n' > "${phylo_out}"/DS_NC2_filarid_trim2.aln
# cat "${phylo_out}"/DS_NC2_non-filarid_trim.aln | tr '\t' '\n' > "${phylo_out}"/DS_NC2_non-filarid_trim2.aln
# "${trimal_cmd}" -in "${phylo_out}"/DS_NC2_non-filarid_trim2.aln -out "${phylo_out}"/DS_NC2_non-filarid_filter_trim.aln -resoverlap 0.70 -seqoverlap 70
# cat "${phylo_out}"/DS_NC2_filarid_trim2.aln "${phylo_out}"/DS_NC2_non-filarid_filter_trim.aln > "${phylo_out}"/DS_NC2_trim_filter.aln
## Change to single-line FASTA
# awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < "${phylo_out}"/DS_NC2_trim_filter.aln > "${phylo_out}"/DS_NC2_trim_filter-single.aln
# mv "${phylo_out}"/DS_NC2_trim_filter-single.aln "${phylo_out}"/DS_CHEMO.aln
### Get IDs and compare lists
# awk 'NR%2==1' "${phylo_out}"/DS_NC2_trim.aln > "${phylo_out}"/trimmed_ids.txt
# awk 'NR%2==1' "${phylo_out}"/DS_CHEMO.aln > "${phylo_out}"/filtered_ids.txt
# grep -v -f "${phylo_out}"/filtered_ids.txt "${phylo_out}"/trimmed_ids.txt > "${phylo_out}"/missing_ids.txt


### MrBayes
# mpirun -np 4 ~/install/MrBayes/src/mb ${local_dir}/NChR/phylo/DS_NC.nxs





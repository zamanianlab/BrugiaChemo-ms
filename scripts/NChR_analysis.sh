#!bin/bash

### Define project directories
boxdr=~/Box\ Sync
proj="Local_50HGI"

gh_dir="${boxdr}/GitHub/${proj}"
local_dir="${boxdr}/GHdata/${proj}"

# Species lists based on genome quality
species_gold="${gh_dir}/auxillary/species_gold.txt"
species_ngf="${gh_dir}/auxillary/species_nongold_filarids.txt"
species_ngo="${gh_dir}/auxillary/species_nongold_other.txt"

gold_dir="${local_dir}/g_genomes"
ngf_dir="${local_dir}/ngf_genomes"
ngo_dir="${local_dir}/ngo_genomes"

### Define script output directories
mkdir "${local_dir}/NChR/"
mkdir "${local_dir}/NChR/g_genomes"
gold_out="${local_dir}/NChR/g_genomes"
mkdir "${local_dir}/NChR/ngf_genomes"
ngf_out="${local_dir}/NChR/ngf_genomes"
mkdir "${local_dir}/NChR/ngo_genomes"
ngo_out="${local_dir}/NChR/ngo_genomes"

##Auxillary Scripts
# Extracting sequences provided list of sequence names and fasta file
seqextract_py="${gh_dir}"/scripts/auxillary/seq_extract.py
# HMMTOP parsing script (filter based on TM range and produce sequences based on TM domains)
HMMTOP_py="${gh_dir}"/scripts/auxillary/HMMTOP_extract.py

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



### Mine Gold Genomes for nematode chemo Rs
### line = species name, iterate through gold genome species names
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


# ##Parse hmm outputs to filter out those where first hit is not NChR hmm, extract sequences of surviving hits
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

# # ##Reciprocal HMMSEARCH of extracted sequences against pfam-a
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		#HMMSEARCH all putative chemosensory genes against db of PFAM hmms
# 		hmmsearch --tblout "${gold_out}"/${line}_rHMM.out --noali "${pfam_HMM}" "${gold_out}"/${line}_NC.fa
# 	done;
# done <"$species_gold"
		

# ###Parse hmm outputs to remove sequences where first hit is not NChR, get list of surviving unique seq ids, extract sequences 
# while IFS= read -r line; do
# 	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
# 		cat "${gold_out}"/${line}_rHMM.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7TM_GPCR_S|srg|sre|srxa/' | sort -k4 -g > "${gold_out}"/${line}_NChitsf.txt
# 		cat "${gold_out}"/${line}_rHMM.out | awk '{print $1 " " $3  " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7TM_GPCR_S|srg|sre|srxa/' | sort -k4 -g  | awk '{print $1}' > "${gold_out}"/${line}_NChitsf_ids.txt
# 		curr_dir=$(dirname "${f}")
#  		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${gold_out}"/${line}_NChitsf_ids.txt "${curr_dir}"/protein.tmp.fa "${gold_out}"/${line}_NCf.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species_gold"

### Remove C. elegans pseudogenes
# rm "${gold_out}"/caenorhabditis_elegans_NCf.fa

## extract pseudogenes from C. elegans GTF
# gtf_parser="${gh_dir}"/scripts/auxillary/gtf_parse.py
# gzcat "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS9.canonical_geneset.gtf.gz > "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS9.canonical_geneset.gtf
# python "${gtf_parser}" "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS9.canonical_geneset.gtf > "${local_dir}"/auxillary/GTF_pseudogenes.txt
# awk '{print $6}' "${local_dir}"/auxillary/GTF_pseudogenes.txt > "${local_dir}"/auxillary/GTF_pseudogenes_Tid.txt
# rm "${gold_dir}"/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS9.canonical_geneset.gtf
## compare list of pseudogenes to list of NCf hits; remove pseudogenes from list of NCf
# grep -Fxv -f "${local_dir}"/auxillary/GTF_pseudogenes_Tid.txt "${gold_out}"/caenorhabditis_elegans_NChitsf_ids.txt > "${gold_out}"/caenorhabditis_elegans_NChitsff_ids.txt


# ## CREATE FILARID GENOME TXT FILE AND RUN ONLY ON THAT
# ###Parse hmm outputs for each species to get list of unique seq ids for R-A, R-P, G, F, A/S (can also use  if ($4 <= 1e-50) in awk)
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

## SUPPLEMENT with curated C. elegans and outgroup GRAFS representatives (drosophila, human) GPCRDB?


# Compare pseudogene list 

######
###### PHYLOGENETIC ANALYSIS
######

phylo_dir="${local_dir}/phylo/NChR"

### Copy sequence files to ../phylo/NemChR directory
while IFS= read -r line; do
	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
		cp "${gold_out}"/${line}_*.fa "${phylo_dir}"
	done;
done <"$species_gold"

# ### Cat and label files for alignment/phylo
cat "${phylo_dir}"/*_NCf.fa > "${phylo_dir}"/All_NCf.fa
cat "${phylo_dir}"/*_R.fa > "${phylo_dir}"/All_R.fa
cat "${phylo_dir}"/*_G.fa > "${phylo_dir}"/All_G.fa
cat "${phylo_dir}"/*_F.fa > "${phylo_dir}"/All_F.fa
cat "${phylo_dir}"/*_AS.fa > "${phylo_dir}"/All_AS.fa

# ### Label each sequence with its predicted family
sed 's/>/>NC:/' "${phylo_dir}"/All_NCf.fa > "${phylo_dir}"/All_NCf_lab.fa
sed 's/>/>R:/' "${phylo_dir}"/All_R.fa > "${phylo_dir}"/All_Rf_lab.fa
sed 's/>/>G:/' "${phylo_dir}"/All_G.fa > "${phylo_dir}"/All_Gf_lab.fa
sed 's/>/>F:/' "${phylo_dir}"/All_F.fa > "${phylo_dir}"/All_Ff_lab.fa
sed 's/>/>AS:/' "${phylo_dir}"/All_AS.fa > "${phylo_dir}"/All_ASf_lab.fa
cat "${phylo_dir}"/All*_lab.fa > "${phylo_dir}"/All.fa

# hmmtop="${gh_dir}"/scripts/auxillary/hmmtop2.1/./hmmtop

# cp "${gh_dir}"/scripts/auxillary/hmmtop_2.1/hmmtop.arch .
# cp "${gh_dir}"/scripts/auxillary/hmmtop_2.1/hmmtop.psv .

outgroup_fa="${gh_dir}/auxillary/outgroup.fa"

# # ### HMMTOP
cat "${phylo_dir}"/All_NCf.fa "${outgroup_fa}" > "${phylo_dir}"/All_NCf_outgroup.fa
# "${hmmtop}" -if="${phylo_dir}"/All_NCf_outgroup.fa -of="${phylo_dir}"/All_NCf_outgroup_hmmtop_output.txt -sf=FAS
# # Parse HHMTOP output to get list of seq ids with >= 5 TM domains <= 10 TM domains
# python "${HMMTOP_py}" "${phylo_dir}"/All_NCf_outgroup_hmmtop_output.txt "${phylo_dir}"/All_NCf_outgroup.fa "${phylo_dir}"/All_NCf_outgroup_TMfiltered.fa

### Align files
# mafft --op 2 --ep 1 --thread 2 --maxiterate 1 "${phylo_dir}"/All_NCf_outgroup_TMfiltered.fa > "${phylo_dir}"/All_NCf_outgroup_TMfiltered_align.fa
mafft --op 2 --ep 1 --thread 2 --maxiterate 1 "${phylo_dir}"/All_NCf_outgroup.fa > "${phylo_dir}"/All_NCf_outgroup_align.fa

# #Trim the alignment
# sed 's/:/_/' ${phylo_dir}/All_TMfiltered_align.fa > ${phylo_dir}/All_TMfiltered_align_rename.fa
# cd /Users/mzamanian/Bioinformatics/Noisy-1.5.12 
# noisy  -v ${phylo_dir}/All_TMfiltered_align_rename.fas

# trimal="${gh_dir}"/scripts/auxillary/trimal/source/./trimal

# "${trimal}" -in "${phylo_dir}"/All_NCf_outgroup_TMfiltered_align.fa -out "${phylo_dir}"/All_NCf_outgroup_TMfiltered_align_trim.fa -gt 0.7

# raxml="${gh_dir}"/scripts/auxillary/raxmlHPC-PTHREADS-SSE3

# "${raxml}" -T 2 -f a -x 12345 -p 12345 -# 100 -m PROTCATAUTO -s "${phylo_dir}"/All_NCf_outgroup_TMfiltered_align_trim.fa -n All_NCf_outgroup.ml.tre
# mv All_NCf_outgroup.ml.tre "${phylo_dir}"/
























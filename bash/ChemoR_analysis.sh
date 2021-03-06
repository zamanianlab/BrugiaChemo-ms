#!/bin/bash

### Preparation
proj="BrugiaChemo-ms"

gh_dir="${GIT_PATH}"/"${proj}"
local_dir="${GIT_DATA}"/"${proj}"

## Species list
species="${gh_dir}"/auxillary/species.txt
# chemoR_families="${gh_dir}"/auxillary/chemoR_families.txt

genome_dir="${local_dir}"/genomes

## Define script output directories
# mkdir "${local_dir}"/ChemoR/
# mkdir "${local_dir}"/ChemoR/output
# mkdir "${local_dir}"/ChemoR/phylo

out="${local_dir}"/ChemoR/output
phylo_out="${local_dir}"/ChemoR/phylo
# mcl_out="${local_dir}"/ChemoR/mcl

## Auxillary Scripts
## extract sequences provided list of sequence names and fasta file
seqextract_py="${gh_dir}"/scripts/auxillary/seq_extract.py
## HMMTOP parsing script (filter based on TM range and produce sequences based on TM domains)
HMMTOP_py="${gh_dir}"/scripts/auxillary/HMMTOP_extract.py
HMMTOP_strict_py="${gh_dir}"/scripts/auxillary/HMMTOP_extract_strict.py
## add species labels to FASTA IDs
change_ID_py="${gh_dir}"/scripts/auxillary/id_change.py
## query WormBase ParaSite API to get paralogues
query_api="${gh_dir}"/scripts/auxillary/json_parser.py
## misc
linearize="${gh_dir}"/scripts/auxillary/linearizefasta.awk
trimal_cmd="${gh_dir}"/scripts/auxillary/trimal/source/./trimal
mafft_cmd="${gh_dir}"/scripts/auxillary/mafft
# muscle_cmd="${gh_dir}"/scripts/auxillary/muscle
# makeblastdb_cmd="${gh_dir}"/scripts/auxillary/makeblastdb
# blastp_cmd="${gh_dir}"/scripts/auxillary/blastp

### Build HMMs

## ChemoRs
# cat "${gh_dir}"/auxillary/pfam_HMMs/GPCR/ChemoR/*hmm > "${gh_dir}"/auxillary/pfam_HMMs/GPCR/ChemoR/ChemoR.hmm
# hmmpress "${gh_dir}"/auxillary/pfam_HMMs/GPCR/ChemoR/ChemoR.hmm
# ChemoR_HMM="${gh_dir}"/auxillary/pfam_HMMs/GPCR/ChemoR/ChemoR.hmm

## GRAFS+
# cat "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/*hmm > "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm
# hmmpress "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm
# GRAFS_HMM="${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm

## Combined (GRAFS+ / ChemoR)
# cat "${gh_dir}"/auxillary/pfam_HMMs/GPCR/Primary/GPCRfams.hmm "${gh_dir}"/auxillary/pfam_HMMs/GPCR/ChemoR/ChemoR.hmm > "${gh_dir}"/auxillary/pfam_HMMs/GPCR/GRAFS_NemChR.hmm
# hmmpress "${gh_dir}"/auxillary/pfam_HMMs/GPCR/GRAFS_NemChR.hmm
# GRAFS_ChemoR_HMM="${gh_dir}"/auxillary/pfam_HMMs/GPCR/GRAFS_ChemoR.hmm

## Prepare Pfam-A HMM db
# mkdir "$local_dir/auxillary
# mkdir "$local_dir/auxillary/HMMs
# wget -nc -O "$local_dir/auxillary/HMMs/Pfam-A.hmm.gz" ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.g
# gzcat "$local_dir/auxillary/HMMs/Pfam-A.hmm.gz" > "$local_dir/auxillary/HMMs/Pfam-A.hmm
# hmmpress "$local_dir/auxillary/HMMs/Pfam-A.hmm
# pfam_HMM="${local_dir}"/auxillary/HMMs/Pfam-A.hmm


############################################################################
###########																											 ###########
###########																											 ###########
###########										ChemoR Identification							 ###########
###########																											 ###########
###########																											 ###########
############################################################################


### Mine for nematode chemo Rs
# while IFS= read -r line; do
# 	for f in "${genome_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 	  	curr_dir=$(dirname "${f}")
# 	  	# echo "${curr_dir}"
# 		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
#
# 		#HMMSEARCH all proteomes against db of All GPCR hmms
# 		hmmsearch --tblout "${out}"/"${line}"_1.out --noali "${GRAFS_ChemoR_HMM}" "${curr_dir}"/protein.tmp.fa
#
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species"


### Parse hmm outputs to filter out those where first hit is not ChemoR or 7tm_1 HMM, extract sequences of surviving hits
# while IFS= read -r line; do
#   for f in "${genome_dir}"/"${line}"/**/*.protein.fa.gz ; do
#     cat "${out}"/"${line}"_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g > "${out}"/"${line}"_1.txt
#     cat "${out}"/"${line}"_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g | awk '{print $1}' > "${out}"/"${line}"_1_ids.txt
#     #Extract these sequences
#     curr_dir=$(dirname "${f}")
#     #echo "${curr_dir}"
#     gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
#     python "${seqextract_py}" "${out}"/"${line}"_1_ids.txt "${curr_dir}"/protein.tmp.fa "${out}"/"${line}"_1.fa
#     rm "${curr_dir}"/protein.tmp.fa
#   done;
# done <"$species"


### Reciprocal HMMSEARCH of extracted sequences against pfam-a
# while IFS= read -r line; do
# 	for f in "${genome_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		#HMMSEARCH all putative chemosensory genes against db of PFAM hmms
# 		hmmsearch --tblout "${out}"/"${line}"_2.out --noali "${pfam_HMM}" "${out}"/"${line}"_1.fa
# 	done;
# done <"$species"


### Parse hmm outputs to remove sequences where first hit is not ChemoR or 7tm_1 HMM, get list of surviving unique IDs, extract sequences
# while IFS= read -r line; do
#   for f in "${genome_dir}"/"${line}"/**/*.protein.fa.gz ; do
#     cat "${out}"/"${line}"_2.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1|7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g > "${out}"/"${line}"_2.txt
#     cat "${out}"/"${line}"_2.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1|7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g  | awk '{print $1}' > "${out}"/"${line}"_2_ids.txt
#     curr_dir=$(dirname "${f}")
#     gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
#     python "${seqextract_py}" "${out}"/"${line}"_2_ids.txt "${curr_dir}"/protein.tmp.fa "${out}"/"${line}"_2.fa
#     rm "${curr_dir}"/protein.tmp.fa
#     # compare the list from _1 to the list from _2 and write out the IDs that were removed
#     grep -v -f "${out}"/"${line}"_2_ids.txt "${out}"/"${line}"_1_ids.txt > "${out}"/"${line}"_2_filtered.txt
#   done;
# done <"$species"


### Parse hmm outputs for each species to get list of unique seq ids for R-A, R-P, G, F, A/S (can also use  if ($4 <= 1e-50) in awk)
# while IFS= read -r line; do
# 	for f in "${genome_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		#Rhodopsin
# 		cat "${out}"/"${line}"_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g > "${out}"/"${line}"_R.txt
# 		cat "${out}"/"${line}"_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g  | awk '{if ($4 <= 1e-0) print $1}' > "${out}"/"${line}"_R_ids.txt
# 		curr_dir=$(dirname "${f}")
#  		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa
# 		python "${seqextract_py}" "${out}"/"${line}"_R_ids.txt "${curr_dir}"/protein.tmp.fa "${out}"/"${line}"_R.fa
# 		#Glutamate
# 		cat "${out}"/"${line}"_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_3/' | sort -k4 -g > "${out}"/"${line}"_G.txt
# 		cat "${out}"/"${line}"_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_3/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${out}"/"${line}"_G_ids.txt
# 		python "${seqextract_py}" "${out}"/"${line}"_G_ids.txt "${curr_dir}"/protein.tmp.fa "${out}"/"${line}"_G.fa
# 		#Adhesion/Secretin
# 		cat "${out}"/"${line}"_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_2/' | sort -k4 -g > "${out}"/"${line}"_AS.txt
# 		cat "${out}"/"${line}"_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_2/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${out}"/"${line}"_AS_ids.txt
# 		python "${seqextract_py}" "${out}"/"${line}"_AS_ids.txt "${curr_dir}"/protein.tmp.fa "${out}"/"${line}"_AS.fa
# 		#Frizzled
# 		cat "${out}"/"${line}"_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/Frizzled/' | sort -k4 -g > "${out}"/"${line}"_F.txt
# 		cat "${out}"/"${line}"_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/Frizzled/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $1}' > "${out}"/"${line}"_F_ids.txt
# 		python "${seqextract_py}" "${out}"/"${line}"_F_ids.txt "${curr_dir}"/protein.tmp.fa "${out}"/"${line}"_F.fa
# 		rm "${curr_dir}"/protein.tmp.fa
# 	done;
# done <"$species"


### Reciprocal blastp of extracted sequences against C. elegans
# while IFS= read -r line; do
#   for f in "${genome_dir}"/"${line}"/**/*.protein.fa.gz ; do
#     # blast filtered ChemoRs against C. elegans proteome, using E-value cutoff
#     cd "${genome_dir}"/caenorhabditis_elegans/PRJNA13758/
#     blastp -query "${out}"/"${line}"_2.fa -db caenorhabditis_elegans.PRJNA13758.WBPS13.protein.fa -out "${out}"/"${line}"_3.out -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
#   done;
# done <"$species"


### Run json_parser.py to populate list of C. elegans ChemoR paralogs to use a comparison for the reciprocal BLAST output
# python "${query_api}"
## use WormBase SimpleMine to convert Gene_IDs to Sequence_IDs (becauase BLASTp uses Sequence_IDs)
## upload list of Gene_IDs that result from the above command
## choose "allow duplicate genes" and "Transcript"
## in Sublime, remove transcript numbers (e.g. N.1) and put each isoform on a new line
## 2177 transcripts


### Remove hits that aren't most similar to a C. elegans ChemoR, extract sequences of surviving hits
# while IFS= read -r line; do
#   for f in "${genome_dir}"/"${line}"/**/*.protein.fa.gz ; do
#     ## Find genes that had no hit at all
#     cat "${out}"/"${line}"_3.out | awk '{print $1}' | grep -wFv -f - "${out}"/"${line}"_2_ids.txt > "${out}"/"${line}"_3_nohits.txt
#     ## Find genes that had a significant hit to a C. elegans ChemoR
#     cat "${out}"/"${line}"_3.out | awk '{print $1 " " $2 " " $7}' | sort -k1,1 -k3,3g | sort -uk1,1 | grep -wF -f "${gh_dir}"/auxillary/ChemoR/simplemine_results.txt | sort -k3 -g > "${out}"/"${line}"_3.txt
#     cat "${out}"/"${line}"_3.txt | awk '{print $1}' | cat - "${out}"/"${line}"_3_nohits.txt > "${out}"/"${line}"_3_ids.txt
#     ## Extract these sequences
#     curr_dir=$(dirname "${f}")
#     echo "${curr_dir}"
#     gzcat -d -k "${f}" > "${curr_dir}"/protein.tmp.fa
#     python "${seqextract_py}" "${out}"/"${line}"_3_ids.txt "${curr_dir}"/protein.tmp.fa "${out}"/"${line}"_3.fa
#     rm "${curr_dir}"/protein.tmp.fa
#     grep -v -f "${out}"/"${line}"_3_ids.txt "${out}"/"${line}"_2_ids.txt > "${out}"/"${line}"_3_filtered.txt
#   done;
# done <"$species"

############################################################################
###########																											 ###########
###########																											 ###########
###########												Phylogenetics									 ###########
###########																											 ###########
###########																											 ###########
############################################################################

### Copy sequence files to ChemoR/phylo directory
# while IFS= read -r line; do
#   for f in "${genome_dir}"/"${line}"/**/*.protein.fa.gz ; do
#     cp "${out}"/"${line}"_3.fa "${phylo_out}"/1/
#   done;
# done <"$species"

## Label each sequence with its species name
# for f in "${phylo_out}"/1/*_3.fa ; do
#   python "${change_ID_py}" "$f";
# done

### Choose one or more representatives from each clade (5354 sequences)
# cat "${phylo_out}"/1/trichinella_spiralis_3_label.fa \
  #   "${phylo_out}"/1/romanomermis_culicivorax_3_label.fa \
  #   "${phylo_out}"/1/syphacia_muris_3_label.fa \
  #   "${phylo_out}"/1/ascaris_suum_3_label.fa \
  #   "${phylo_out}"/1/toxocara_canis_3_label.fa  \
  #   "${phylo_out}"/1/brugia_malayi_3_label.fa \
  #   "${phylo_out}"/1/onchocerca_volvulus_3_label.fa \
  #   "${phylo_out}"/1/strongyloides_ratti_3_label.fa \
  #   "${phylo_out}"/1/rhabditophanes_kr3021_3_label.fa \
  #   "${phylo_out}"/1/meloidogyne_hapla_3_label.fa \
  #   "${phylo_out}"/1/panagrellus_redivivus_3_label.fa \
  #   "${phylo_out}"/1/haemonchus_contortus_3_label.fa \
  #   "${phylo_out}"/1/nippostrongylus_brasiliensis_3_label.fa \
  #   "${phylo_out}"/1/angiostrongylus_cantonensis_3_label.fa \
  #   "${phylo_out}"/1/dictyocaulus_viviparus_3_label.fa \
  #   "${phylo_out}"/1/necator_americanus_3_label.fa \
  #   "${phylo_out}"/1/ancylostoma_caninum_3_label.fa \
  #   "${phylo_out}"/1/pristionchus_pacificus_3_label.fa > \
  # "${phylo_out}"/2/down_sampled_1.fa

## HMMTOP
# cd "${gh_dir}"/scripts/auxillary/hmmtop_2.1/
# ./hmmtop -if="${phylo_out}"/2/down_sampled_1.fa -of="${phylo_out}"/2/down_sampled_hmmtop_output.txt -sf=FAS

### Parse HHMTOP output to get FASTA file of non-filarid sequences with exactly 7 TMs; extract entire sequence (1496 with 7 TMs)
# python "${HMMTOP_strict_py}" "${phylo_out}"/2/down_sampled_hmmtop_output.txt "${phylo_out}"/2/down_sampled_1.fa "${phylo_out}"/2/down_sampled_2.fa

### prepare C. elegans ChemoR fasta file
# python "${seqextract_py}" "${gh_dir}"/auxillary/ChemoR/simplemine_results.txt "${genome_dir}"/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS13.protein.fa "${gh_dir}"/auxillary/ChemoR/simplemine_results.fa
# awk -f "${linearize}" < "${gh_dir}"/auxillary/ChemoR/simplemine_results.fa > "${gh_dir}"/auxillary/ChemoR/simplemine_results_linear.fa
# sed -i -E 's/\t/\n/g' "${gh_dir}"/auxillary/ChemoR/simplemine_results_linear.fa

### TM prediction
# cd "${gh_dir}"/scripts/auxillary/hmmtop_2.1/
# ./hmmtop -if="${gh_dir}"/auxillary/ChemoR/simplemine_results_linear.fa -of="${gh_dir}"/auxillary/ChemoR/simplemine_results_hmmtop.txt -sf=FAS
# python "${HMMTOP_strict_py}" "${gh_dir}"/auxillary/ChemoR/simplemine_results_hmmtop.txt "${gh_dir}"/auxillary/ChemoR/simplemine_results_linear.fa "${gh_dir}"/auxillary/ChemoR/simplemine_results_3.fa
# awk -f "${linearize}" < "${gh_dir}"/auxillary/ChemoR/simplemine_results_3.fa > "${gh_dir}"/auxillary/ChemoR/simplemine_results_4.fa
# gsed -i -E 's/\t/\n/g' "${gh_dir}"/auxillary/ChemoR/simplemine_results_4.fa
# gsed -i -E 's/^>/>cele-/g' "${gh_dir}"/auxillary/ChemoR/simplemine_results_4.fa
# mv "${gh_dir}"/auxillary/ChemoR/simplemine_results_4.fa "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa


### Create a separate fasta file for each family
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'sra' > "${phylo_out}"/3/sra.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srab' > "${phylo_out}"/3/srab.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srh' > "${phylo_out}"/3/srh.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'str' > "${phylo_out}"/3/str.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'sri' > "${phylo_out}"/3/sri.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srd' > "${phylo_out}"/3/srd.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srj' > "${phylo_out}"/3/srj.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srm' > "${phylo_out}"/3/srm.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srn' > "${phylo_out}"/3/srn.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'sre' > "${phylo_out}"/3/sre.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srb' > "${phylo_out}"/3/srb.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srx' > "${phylo_out}"/3/srx.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srt' > "${phylo_out}"/3/srt.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srg' > "${phylo_out}"/3/srg.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'sru' > "${phylo_out}"/3/sru.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'sro' > "${phylo_out}"/3/sro.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srv' > "${phylo_out}"/3/srv.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srxa' > "${phylo_out}"/3/srxa.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srw' > "${phylo_out}"/3/srw.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srz' > "${phylo_out}"/3/srz.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srbc' > "${phylo_out}"/3/srbc.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srsx' > "${phylo_out}"/3/srsx.celeg.fa
# cat "${gh_dir}"/auxillary/ChemoR/caenorhabditis_elegans_curated.fa | grep -wF -A1 'srr' > "${phylo_out}"/3/srr.celeg.fa

### align family fasta files
# for f in "${phylo_out}"/3/s*.celeg.fa; do
#   "mafft" --auto "$f" > "$f.aln";
# done

# for f in "${phylo_out}"/3/*.fa.aln ; do
#   mv "$f" "${f/.fa.aln/.aln}";
# done

### align family profiles to each other, sequentially
# muscle -profile -in1 "${phylo_out}"/3/sra.celeg.aln -in2 "${phylo_out}"/3/srab.celeg.aln -out "${phylo_out}"/3/2.aln
# muscle -profile -in1 "${phylo_out}"/3/2.aln -in2 "${phylo_out}"/3/srh.celeg.aln -out "${phylo_out}"/3/3.aln
# muscle -profile -in1 "${phylo_out}"/3/3.aln -in2 "${phylo_out}"/3/str.celeg.aln -out "${phylo_out}"/3/4.aln
# muscle -profile -in1 "${phylo_out}"/3/4.aln -in2 "${phylo_out}"/3/sri.celeg.aln -out "${phylo_out}"/3/5.aln
# muscle -profile -in1 "${phylo_out}"/3/5.aln -in2 "${phylo_out}"/3/srd.celeg.aln -out "${phylo_out}"/3/6.aln
# muscle -profile -in1 "${phylo_out}"/3/6.aln -in2 "${phylo_out}"/3/srj.celeg.aln -out "${phylo_out}"/3/7.aln
# muscle -profile -in1 "${phylo_out}"/3/7.aln -in2 "${phylo_out}"/3/sre.celeg.aln -out "${phylo_out}"/3/8.aln
# muscle -profile -in1 "${phylo_out}"/3/8.aln -in2 "${phylo_out}"/3/srb.celeg.aln -out "${phylo_out}"/3/9.aln
# muscle -profile -in1 "${phylo_out}"/3/9.aln -in2 "${phylo_out}"/3/srx.celeg.aln -out "${phylo_out}"/3/10.aln
# muscle -profile -in1 "${phylo_out}"/3/10.aln -in2 "${phylo_out}"/3/srt.celeg.aln -out "${phylo_out}"/3/11.aln
# muscle -profile -in1 "${phylo_out}"/3/11.aln -in2 "${phylo_out}"/3/srg.celeg.aln -out "${phylo_out}"/3/12.aln
# muscle -profile -in1 "${phylo_out}"/3/12.aln -in2 "${phylo_out}"/3/sru.celeg.aln -out "${phylo_out}"/3/13.aln
# muscle -profile -in1 "${phylo_out}"/3/13.aln -in2 "${phylo_out}"/3/srxa.celeg.aln -out "${phylo_out}"/3/14.aln
# muscle -profile -in1 "${phylo_out}"/3/14.aln -in2 "${phylo_out}"/3/srw.celeg.aln -out "${phylo_out}"/3/15.aln
# muscle -profile -in1 "${phylo_out}"/3/15.aln -in2 "${phylo_out}"/3/srz.celeg.aln -out "${phylo_out}"/3/16.aln
# muscle -profile -in1 "${phylo_out}"/3/16.aln -in2 "${phylo_out}"/3/srbc.celeg.aln -out "${phylo_out}"/3/17.aln
# muscle -profile -in1 "${phylo_out}"/3/17.aln -in2 "${phylo_out}"/3/srsx.celeg.aln -out "${phylo_out}"/3/18.aln
# muscle -profile -in1 "${phylo_out}"/3/18.aln -in2 "${phylo_out}"/3/srv.celeg.aln -out "${phylo_out}"/3/19.aln
# muscle -profile -in1 "${phylo_out}"/3/19.aln -in2 "${phylo_out}"/3/srm.celeg.aln -out "${phylo_out}"/3/20.aln
# muscle -profile -in1 "${phylo_out}"/3/20.aln -in2 "${phylo_out}"/3/srn.celeg.aln -out "${phylo_out}"/3/21.aln
# muscle -profile -in1 "${phylo_out}"/3/21.aln -in2 "${phylo_out}"/3/sro.celeg.aln -out "${phylo_out}"/3/22.aln
# muscle -profile -in1 "${phylo_out}"/3/22.aln -in2 "${phylo_out}"/3/srr.celeg.aln -out "${phylo_out}"/3/23.aln

### add non-C.elegans representatives to the alignment
# mafft --reorder --thread 8 --addfull "${phylo_out}"/2/down_sampled_2.fa --keeplength "${phylo_out}"/3/23.aln > "${phylo_out}"/4/ds_reps.aln
### trim and filter
# "${trimal_cmd}" -gt 0.7 -in "${phylo_out}"/4/ds_reps.aln -out "${phylo_out}"/4/ds_reps_2.aln
# "${trimal_cmd}" -resoverlap 0.70 -seqoverlap 70 -in "${phylo_out}"/4/ds_reps_2.aln -out "${phylo_out}"/4/ds_reps_3.aln
# ### Change to single-line FASTA
# awk -f "${linearize}" < "${phylo_out}"/4/ds_reps_2.aln > "${phylo_out}"/4/ds_reps_2_linear.aln
# awk -f "${linearize}" < "${phylo_out}"/4/ds_reps_3.aln | sed '/^$/d' > "${phylo_out}"/4/ds_reps_3_linear.aln
### Get IDs and compare lists
# cat "${phylo_out}"/4/ds_reps_2_linear.aln | tr '\t' '\n' | awk 'NR%2==1' > "${phylo_out}"/4/ds_reps_2_ids.txt
# cat "${phylo_out}"/4/ds_reps_3_linear.aln | tr '\t' '\n' | awk 'NR%2==1' > "${phylo_out}"/4/ds_reps_3_ids.txt
# grep -v -f "${phylo_out}"/4/ds_reps_3_ids.txt "${phylo_out}"/4/ds_reps_2_ids.txt > "${phylo_out}"/4/ds_reps_3_filtered.txt

### ML tree on server
# iqtree -s ../4/ds_reps_3.aln -nt 4 -alrt 1000 -bb 1000

################################################################################################################################
###########																											 ###########
###########																											 ###########
###########												Heatmap														 ###########
###########																											 ###########
###########																											 ###########
################################################################################################################################


# find -name '*.protein.gz' -exec gzip -d {}
# find -name '*protein.fa' -type f | xargs cat > all_protein.fa
# cat "${phylo_out}"/6/tree_clades.csv | awk -F , '{print $3}' | sed 's/"//g' | tail -n +2 > "${phylo_out}"/6/tree_clades_ids.txt
# python "${seqextract_py}" "${phylo_out}"/6/tree_clades_ids.txt "${local_dir}"/genomes/all_protein.fa "${phylo_out}"/6/tree_clades.fa
# makeblastdb  -dbtype prot -in "${phylo_out}"/6/tree_clades.fa -out "${phylo_out}"/6/tree_clades.db

### BLAST species not included in original tree against families
# while IFS= read -r line; do
# 	for f in "${genome_dir}"/"${line}"/**/*.protein.fa.gz ; do
# 		cd "${phylo_out}"/6/
# 		blastp -query "${phylo_out}"/1/"${line}"_3_label.fa -db tree_clades.db -out "${phylo_out}"/6/"${line}".blastout -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# 	done;
# done <"$species"

### filter BLAST results
# while IFS= read -r line; do
# 	for f in "${genome_dir}"/"${line}"/**/*.protein.fa.gz ; do
#     cat "${phylo_out}"/6/"${line}".blastout | sort -k1,1 -k7,7g | sort -uk1,1 > "${phylo_out}"/6/"${line}"_2.blastout
#   done;
# done <"$species"

## concatenate blast results and remove self to self HSPs
# cat "${phylo_out}"/6/*_2.blastout | sed 's/-/\t/' | awk '$2!=$3' > "${phylo_out}"/6/ChemoR.blastout

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

## get sequences for families, align and trim
# while IFS= read -r line; do
# 	python "${seqextract_py}" "${phylo_out}"/ML/families/"${line}"_ids.txt "${local_dir}"/a_genomes/all.fa "${phylo_out}"/ML/families/"${line}".fa
# 	"${mafft_cmd}" --thread 4 --reorder --auto "${phylo_out}"/ML/families/"${line}".fa > "${phylo_out}"/ML/families/"${line}".aln
# 	"${trimal_cmd}" -gt 0.7 -in "${phylo_out}"/ML/families/"${line}".aln -out "${phylo_out}"/ML/families/"${line}".trim.aln
# 	"${trimal_cmd}" -resoverlap 0.70 -seqoverlap 70 -in "${phylo_out}"/ML/families/"${line}".trim.aln -out "${phylo_out}"/ML/families/"${line}".trim.filter.aln
# done <"$chemoR_families"

## align family profiles to each other to make superfamily alignments
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

# while IFS= read -r line; do
# 	cd "${phylo_out}"/ML/families/
# 	/home/BIOTECH/zamanian/install/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 4 -f a -x 12345 -p 12345 -# 100 -m PROTGAMMAVT -s "${line}.trim.filter.aln" -n "${line}".tre;
# done <"$chemoR_families"

## find sequences that were filtered out
# while IFS= read -r line; do
# 	perl -pe '/^>/ ? print "\n" : chomp' "${phylo_out}"/ML/families/"${line}".trim.aln | tail -n +2 > "${phylo_out}"/ML/families/"${line}".trim-single.aln
# 	perl -pe '/^>/ ? print "\n" : chomp' "${phylo_out}"/ML/families/"${line}".trim.filter.aln | tail -n +2 > "${phylo_out}"/ML/families//"${line}".trim.filter-single.aln
# 	cat "${phylo_out}"/ML/families/"${line}".trim-single.aln | tr '\t' '\n' | awk 'NR%2==1' > "${phylo_out}"/ML/families/"${line}".trimmed_ids.txt
# 	cat "${phylo_out}"/ML/families/"${line}".trim.filter-single.aln | tr '\t' '\n' | awk 'NR%2==1' > "${phylo_out}"/ML/families/"${line}".filtered_ids.txt
# 	grep -v -f "${phylo_out}"/ML/families/"${line}".filtered_ids.txt "${phylo_out}"/ML/families/"${line}".trimmed_ids.txt > "${phylo_out}"/ML/families/"${line}".missing_ids.txt
# done <"$chemoR_families"

# rm "${phylo_out}"/ML/families/*single*

# make a blast db using the sequences that weren't filtered above
# do that by getting the headers, removing the preceding species information, adding a newline to the end of the file, getting the full-length sequence, and making a db
# while IFS= read -r line; do
# 	cd "${phylo_out}"/ML/families
# 	cat "${phylo_out}"/ML/families/"${line}".trim.filter.aln | grep -P '>.*' | sed 's/>//' | sed '$a\' | python "${seqextract_py}" /dev/stdin "${local_dir}"/a_genomes/all.fa "${phylo_out}"/ML/families/"${line}".db.fa
# 	"${makeblastdb_cmd}" -dbtype prot -in "${line}".db.fa -out "${line}".db
# done <"$chemoR_families"

## for those that were removed, blast against families to see where they *would* be placed in tree, if they weren't filtered
# while IFS= read -r line; do
# 	cat "${phylo_out}"/ML/families/"${line}".missing_ids.txt | sed 's/>//' | python "${seqextract_py}" /dev/stdin "${local_dir}"/a_genomes/all.fa "${phylo_out}"/ML/families/"${line}".missing.fa
# 	cd "${phylo_out}"/ML/families/
# 	"${blastp_cmd}" -query "${line}".missing.fa -db "${line}".db -out "${line}".missing.blastout -num_threads 4 -evalue 0.01 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
# done <"$chemoR_families"

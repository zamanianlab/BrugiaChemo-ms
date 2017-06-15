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

## NChR
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

### Mine Gold Genomes for nematode chemo Rs
#line = species name, iterate through gold genome species names
while IFS= read -r line; do
	for f in "${gold_dir}"/${line}/**/*.protein.fa.gz ; do
	  	curr_dir=$(dirname "${f}")
	  	echo "${curr_dir}"
		gzcat "${f}" > "${curr_dir}"/protein.tmp.fa

		#HMMSCAN all proteomes against db of All GPCR hmms
		hmmscan --tblout "${gold_out}"/${line}_hits.out --noali "${GRAFS_NemChR_HMM}" "${curr_dir}"/protein.tmp.fa 

		rm "${curr_dir}"/protein.tmp.fa
	done;
done <"$species_gold"




		# #Parse hmm outputs for each species to get list of unique seq ids for NemChRs (NC)
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '!/Frizzled|7tm_1|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g > ${out_dir}/${line}_NChits.txt
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '!/Frizzled|7tm_1|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g  | awk '{print $3}' > ${out_dir}/${line}_NChits_ids.txt
		# #Extract these sequences
		# python ${seqextract_py} ${out_dir}/${line}_NChits_ids.txt ${curr_dir}/protein.tmp.fa ${out_dir}/${line}_NC.fa

		# #Reciprocal HMMSCAN of extracted sequences against pfam-a
		# hmmscan --tblout ${out_dir}/${line}_rHMM.out --noali ${pfam_hmm} ${out_dir}/${line}_NC.fa

		# #Parse hmm outputs for each species to get list of unique seq ids + extract sequences + TM >=5 filter
		# cat ${out_dir}/${line}_rHMM.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/7TM_GPCR_S|srg|sre|srxa/' | sort -k4 -g > ${out_dir}/${line}_NChitsf.txt
		# cat ${out_dir}/${line}_rHMM.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/7TM_GPCR_S|srg|sre|srxa/' | sort -k4 -g  | awk '{print $3}' > ${out_dir}/${line}_NChitsf_ids.txt
		# python ${seqextract_py} ${out_dir}/${line}_NChitsf_ids.txt ${curr_dir}/protein.tmp.fa ${out_dir}/${line}_NCf.fa

		# #Parse hmm outputs for each species to get list of unique seq ids for R-A, R-P, G, F, A/S (can also use  if ($4 <= 1e-50) in awk)
		# cd /Users/mzamanian/Bioinformatics/hmmtop_2.1/
		# #Rhodopsin
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/7tm_1/' | sort -k4 -g > ${out_dir}/${line}_Rhits.txt
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/7tm_1/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $3}' > ${out_dir}/${line}_Rhits_ids.txt
		# python ${seqextract_py} ${out_dir}/${line}_Rhits_ids.txt ${curr_dir}/protein.tmp.fa ${out_dir}/${line}_R.fa
		# # Glutamate
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/7tm_3/' | sort -k4 -g > ${out_dir}/${line}_Ghits.txt
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/7tm_3/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $3}' > ${out_dir}/${line}_Ghits_ids.txt
		# python ${seqextract_py} ${out_dir}/${line}_Ghits_ids.txt ${curr_dir}/protein.tmp.fa ${out_dir}/${line}_G.fa
		# # Adhesion/Secretin
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/7tm_2/' | sort -k4 -g > ${out_dir}/${line}_AShits.txt
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/7tm_2/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $3}' > ${out_dir}/${line}_AShits_ids.txt
		# python ${seqextract_py} ${out_dir}/${line}_AShits_ids.txt ${curr_dir}/protein.tmp.fa ${out_dir}/${line}_AS.fa
		# # Frizzled
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/Frizzled/' | sort -k4 -g > ${out_dir}/${line}_Fhits.txt
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/Frizzled/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $3}' > ${out_dir}/${line}_Fhits_ids.txt
		# python ${seqextract_py} ${out_dir}/${line}_Fhits_ids.txt ${curr_dir}/protein.tmp.fa ${out_dir}/${line}_F.fa
		# # Insect Odorant/ChemoR
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/Frizzled/' | sort -k4 -g > ${out_dir}/${line}_IOChits.txt
		# cat ${out_dir}/${line}_hits.out | awk '{print $1 " " $2 " " $3 " " $5}' | awk '!/#/' | sort -uk3,3 | awk '/Frizzled/' | sort -k4 -g  | awk '{if ($4 <= 1e-50) print $3}' > ${out_dir}/${line}_IOChits_ids.txt
		# python ${seqextract_py} ${out_dir}/${line}_IOChits_ids.txt ${curr_dir}/protein.tmp.fa ${out_dir}/${line}_IOC.fa



######
###### PHYLOGENETIC ANALYSIS
######

# phylo_dir="/Users/mzamanian/GitHub/ChemoR-sh/phylo/NemChR"

# # # Copy sequence files to ../phylo/NemChR directory
# # while IFS= read -r line; do
# # 	cp ${genome_dir}/${line}/**/NChR/${line}_*.fa ${phylo_dir}
# # done <"$species_gold"

# # # Cat and label files for alignment/phylo
# # cat ${phylo_dir}/*_NCf.fa > ${phylo_dir}/All_NCf.fa
# # cat ${phylo_dir}/c*_elegans_R.fa ${phylo_dir}/a*_suum_R.fa ${phylo_dir}/b*_malayi_R.fa ${phylo_dir}/s*_ratti_R.fa ${phylo_dir}/h*_contortus_R.fa ${phylo_dir}/t*_canis_R.fa > ${phylo_dir}/All_R.fa
# # cat ${phylo_dir}/c*_elegans_F.fa ${phylo_dir}/a*_suum_F.fa ${phylo_dir}/b*_malayi_F.fa ${phylo_dir}/s*_ratti_F.fa ${phylo_dir}/h*_contortus_F.fa ${phylo_dir}/t*_canis_F.fa > ${phylo_dir}/All_F.fa
# # cat ${phylo_dir}/c*_elegans_G.fa ${phylo_dir}/a*_suum_G.fa ${phylo_dir}/b*_malayi_G.fa ${phylo_dir}/s*_ratti_G.fa ${phylo_dir}/h*_contortus_G.fa ${phylo_dir}/t*_canis_G.fa > ${phylo_dir}/All_G.fa
# # cat ${phylo_dir}/c*_elegans_AS.fa ${phylo_dir}/a*_suum_AS.fa ${phylo_dir}/b*_malayi_AS.fa ${phylo_dir}/s*_ratti_AS.fa ${phylo_dir}/h*_contortus_AS.fa ${phylo_dir}/t*_canis_AS.fa > ${phylo_dir}/All_AS.fa

# # # Label each sequence with its predicted family
# # sed 's/>/>NC:/' ${phylo_dir}/All_NCf.fa > ${phylo_dir}/All_NCf_lab.fa
# # sed 's/>/>R:/' ${phylo_dir}/All_R.fa > ${phylo_dir}/All_Rf_lab.fa
# # sed 's/>/>G:/' ${phylo_dir}/All_G.fa > ${phylo_dir}/All_Gf_lab.fa
# # sed 's/>/>F:/' ${phylo_dir}/All_F.fa > ${phylo_dir}/All_Ff_lab.fa
# # sed 's/>/>AS:/' ${phylo_dir}/All_AS.fa > ${phylo_dir}/All_ASf_lab.fa
# # cat ${phylo_dir}/All*_lab.fa > ${phylo_dir}/All.fa

# HMMextract_py="/Users/mzamanian/GitHub/ChemoR-sh/scripts/HMMTOP_extract.py"

# # # HMMTOP 
# # cd /Users/mzamanian/Bioinformatics/hmmtop_2.1/
# # ./hmmtop -if=${phylo_dir}/All.fa -of=${phylo_dir}/hmm_output.txt -sf=FAS
# # Parse HHMTOP output to get list of seq ids with >= 5 TM domains <= 10 TM domains
# #python ${HMMextract_py} ${phylo_dir}/hmm_output.txt ${phylo_dir}/All.fa ${phylo_dir}/All_TMfiltered.fa

# # #Align files
# #mafft --op 2 --ep 1 --thread 2 --maxiterate 1 ${phylo_dir}/All_TMfiltered.fa > ${phylo_dir}/All_TMfiltered_align.fa


# #Trim the alignment
# sed 's/:/_/' ${phylo_dir}/All_TMfiltered_align.fas > ${phylo_dir}/All_TMfiltered_align_rename.fas
# cd /Users/mzamanian/Bioinformatics/Noisy-1.5.12 
# noisy  -v ${phylo_dir}/All_TMfiltered_align_rename.fas

#trimal -in ${phylo_dir}/All_TMfiltered_align_rename -out ${phylo_dir}/All_TMfiltered_align_trim.nex -nexus -gt 0.9 -cons 2




#####

##old approach to extracting hmmtop threshold hits
#cat ${phylo_dir}/hmm_output.txt | awk '{if ($(NF-14) >= 5) print $3}' > ${phylo_dir}/hmm_output_filter.txt
#python ${seqextract_py} ${curr_dir}/${line}_NC2_HMMTOPf.txt ${curr_dir}/protein.tmp.fa ${curr_dir}/${line}_NCf.fa



#Blast additional option
#/usr/local/ncbi/blast/bin/blastp -query ${curr_dir}/${line}_nemchr.fa -db nr -remote -entrez_query "Caenorhabditis elegans[Organism]" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" -max_target_seqs 5 -out ${curr_dir}/${line}_nemchr_blast.txt

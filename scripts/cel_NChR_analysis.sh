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
# mkdir "${local_dir}/NChR/phylo/c_elegans/"

gold_out="${local_dir}/NChR/g_genomes"
ngo_out="${local_dir}/NChR/ngo_genomes"
ngf_out="${local_dir}/NChR/ngf_genomes"
phylo_out="${local_dir}/NChR/phylo"
mcl_out="${local_dir}/NChR/mcl"
cel_out="${local_dir}/NChR/phylo/c_elegans/"

### Start C. elegans NChR re-alignment (following protocol from https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-6-42)
# mv "${phylo_out}"/caenorhabditis_elegans_NCf_label.fa "${cel_out}"/
## Align files
# einsi --thread 8 "${cel_out}"/caenorhabditis_elegans_NCf_label.fa > "${cel_out}"/caenorhabditis_elegans_NCf_label.aln
## Email when complete
# mailx -s "Alignment complete!" njwheeler@wisc.edu <<< "The alignment of caenorhabditis_elegans_NCf_label has successfully completed."
## Trim and filter (filters out 120 sequences with < 70% overlap with the rest of the alignment)
trimal_cmd="${gh_dir}"/scripts/auxillary/trimal/source/./trimal
# "${trimal_cmd}" -in "${cel_out}"/caenorhabditis_elegans_NCf_label.aln -out "${cel_out}"/caenorhabditis_elegans_NCf_trim.aln -gt 0.7 
# "${trimal_cmd}" -in "${cel_out}"/caenorhabditis_elegans_NCf_trim.aln -out "${cel_out}"/caenorhabditis_elegans_NCf_trim_filter.aln -resoverlap 0.70 -seqoverlap 70
## Make single-line FASTA and get IDs to compare lists
# awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < "${cel_out}"/caenorhabditis_elegans_NCf_trim.aln > "${cel_out}"/caenorhabditis_elegans_NCf_trim-single.aln
# awk '/^>/ {printf("%s%s\n",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < "${cel_out}"/caenorhabditis_elegans_NCf_trim_filter.aln > "${cel_out}"/caenorhabditis_elegans_NCf_trim_filter-single.aln
# awk 'NR%2==1' "${cel_out}"/caenorhabditis_elegans_NCf_trim-single.aln > "${cel_out}"/trimmed_ids.txt
# awk 'NR%2==1' $"${cel_out}"/caenorhabditis_elegans_NCf_trim_filter-single.aln > "${cel_out}"/filtered_ids.txt
# grep -v -f "${cel_out}"/filtered_ids.txt "${cel_out}"/trimmed_ids.txt > "${cel_out}"/missing_ids.txt

### ML tree inference
/home/BIOTECH/zamanian/install/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -T 4 -f a -x 12345 -p 12345 -# 100 -m PROTGAMMAAUTO -s "${cel_out}"/caenorhabditis_elegans_NCf_trim_filter.aln -n "${cel_out}"/caenorhabditis_elegans_ML
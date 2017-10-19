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
mv "${phylo_out}"/caenorhabditis_elegans_NCf_label.fa "${cel_out}"/
### Align files
einsi --thread 8 "${cel_out}"/caenorhabditis_elegans_NCf_label.fa > "${cel_out}"/caenorhabditis_elegans_NCf_label.aln
### Email when complete
mailx -s "Alignment complete!" njwheeler@wisc.edu <<< "The alignment of caenorhabditis_elegans_NCf_label has successfully completed."
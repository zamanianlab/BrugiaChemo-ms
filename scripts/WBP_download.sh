#!bin/bash

### Preparation
proj="50HGI"

gh_dir="${GIT_PATH}"/"${proj}"
local_dir="${GIT_DATA}"/"${proj}"

### Define wormbase source links and download genomes

wbp_prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS13/species"

# Species list based on genome quality
species="${gh_dir}/auxillary/species.txt"

# Create folders corresponding to genome destinations
mkdir "${local_dir}"
mkdir "${local_dir}/genomes"

# Added -N to only download newer versions
while IFS= read -r line
do
	species_dl="${wbp_prefix}/$line/"
	printf ${species_dl}"\n"
	wget -nc -r -nH --cut-dirs=7 --no-parent --reject="index.html*" -A 'canonical_geneset.gtf.gz','genomic.fa.gz','protein.fa.gz','CDS_transcripts.fa.gz' $species_dl -P "${local_dir}/genomes"
done <"$species"

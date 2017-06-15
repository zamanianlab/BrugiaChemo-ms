#!bin/bash

### Define project directories
boxdr=~/Box\ Sync
proj="Local_50HGI"

gh_dir="${boxdr}/GitHub/${proj}"
local_dir="${boxdr}/GHdata/${proj}"

### Define wormbase source links and download genomes

wbp_prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS9/species"

# Species lists based on genome quality
species_gold="${gh_dir}/auxillary/species_gold.txt"
species_ngf="${gh_dir}/auxillary/species_nongold_filarids.txt"
species_ngo="${gh_dir}/auxillary/species_nongold_other.txt"

# Create folders corresponding to genome destinations
mkdir "${local_dir}/g_genomes"
mkdir "${local_dir}/ngf_genomes"
mkdir "${local_dir}/ngo_genomes"

# Download gold genomes -> g_genomes folder
# Added -N to only download newer versions
while IFS= read -r line
do
	species_dl="${wbp_prefix}/$line/"
	printf ${species_dl}"\n"
	wget -nc -r -nH --cut-dirs=7 --no-parent --reject="index.html*" -A 'canonical_geneset.gtf.gz','genomic.fa.gz','protein.fa.gz','CDS_transcripts.fa.gz' $species_dl -P "${local_dir}/g_genomes"
done <"$species_gold"

## Download remaining non-gold filarids -> ngf_genomes folder
# Added -N to only download newer versions
while IFS= read -r line
do
	species_dl="${wbp_prefix}/$line/"
	printf ${species_dl}"\n"
	wget -nc -r -nH --cut-dirs=7 --no-parent --reject="index.html*" -A 'canonical_geneset.gtf.gz','genomic.fa.gz','protein.fa.gz','CDS_transcripts.fa.gz' $species_dl -P "${local_dir}/ngf_genomes"
done <"$species_ngf"

# Download remaining non-gold other -> ngo_genomes folder
# Added -N to only download newer versions
while IFS= read -r line
do
	species_dl="${wbp_prefix}/$line/"
	printf ${species_dl}"\n"
	wget -nc -r -nH --cut-dirs=7 --no-parent --reject="index.html*" -A 'canonical_geneset.gtf.gz','genomic.fa.gz','protein.fa.gz','CDS_transcripts.fa.gz' $species_dl -P "${local_dir}/ngo_genomes"
done <"$species_ngo"
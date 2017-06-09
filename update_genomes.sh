#!/bin/bash

# wget command to populate 50HGI data
#cd ../50HGI_data
#wget -r -nH -nd -np ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/*

#lftp
HOST='ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/'
SOURCEFOLDER='species'

lftp -f "
open $HOST
mirror -c -P -v $SOURCEFOLDER ../50HGI_data
bye
"
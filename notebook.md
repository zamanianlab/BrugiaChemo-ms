notebook.md

## Notes

This notebook was began after much of the development of the ChemoR pipeline had been completed and when a significant re-write was undertaken to begin to prepare for publication. The following contains a list of things that needs to be done, and the approaches attempted to achieve these goals.

## Tasks

### Organization

1. Begin to organize with the "\_N" nomenclature to keep output files organized
2. Add better comments in the .sh scripts
3. Rewrite in Nextflow


### Improvements

1. Remove all reliance on the `-max_target_seqs` option in `blastp` commands
2. Improve the final reciprocal BLAST filtering strategy
  - It currently relies upon a `grep` comparison with the output from the script on the *C. elegans* genome
  - It should rely upon a list of hard-coded, curated genes
  - It also removes putative ChemoRs that **do not** have a hit to a *C. elegans* gene, which is a bug
  - It should keep those putative ChemoRs and even focus on them
3. Add a search using putative ChemoRs against the genome that they come from (to get diverged sequences)
4. Be explicit about Gene_ID and Transcript_ID, and make a master data frame with everything
5. TM prediction on everything prior to phylogenetic inference
6. Change what species go into the tree
7. When there are isoforms, keep the longest sequence


### Rolling Notebook

#### 2019-05-02

- began this notebook  
- met with Mostafa to plan ahead for publication/preprint  
- updated `WBP_download.sh` to not discriminate between genome "quality" and include everything in the same `genomes/` directory in `GHdata`
- ran `sh WBP_download.sh` to update the genomes
- deleted `g_genomes`, `ngf_genomes`, `f_genomes`, `ngo_genomes`, `a_genomes`, and `p_genomes`
- keep `7tm_1` and ChemoR hits in the first `hmmsearch` and the reciprocal search
- on line 164 in `ChemoR_analysis.sh`, the bit where the reciprocal BLAST output is compared to *C. elegans* ChemoRs

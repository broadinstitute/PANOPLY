ClumpsPTM requires access to structural information, in order to determine the spatial arrangement of PTMs. Currently, the PANOPLY ClumpsPTM workflow requires access to the Protein Data Bank (PDB); it expects to be provided a Google Bucket address, containing the PDB archive as a set of tar files. The provided script `update_tar_pdb.sh` can be used to generate this. The script first syncs the PDB archive locally, then chunks the directory into a set of tar files, and finally pushes those tarred files to a Google Bucket.

The script also sources two additional references files:
 * The current Uniprot FASTA sequences `uniprot_sprot.fasta` for BLASTing to Uniprot accession IDs
 * The SIFTS database `pdb_chain_uniprot.tsv` for mapping from Uniprot to PDB IDs



- - - ADDITIONAL NOTES - - -

PDBs can be pulled from one of two places, depending on which rsync lines are left uncommented:
 * Current PDB-- Up-to-date PDB archive for most-current structural data (rsync.rcsb.org::ftp_data/structures/divided/pdb/)
 * Frozen PDB-- Snapshot of the PDB frozen to a specific date; allows for reproducible results. The current snapshot is frozen to 2025; it is recommended to update this if you are running in a later year. (snapshots.pdbj.org::20250101/pub/pdb/data/structures/divided/pdb/)

In addition to a full copy of the PDB-- this script additionally creates a subsetted copy of the PDB for testing purposes. In this testing set, most tar-files are empty and only a few contain structural data; this directory is very lightweight and can be helpful for debugging the PANOPLY ClumpsPTM workflow, since the smaller tar-files requires significantly less time to localize.



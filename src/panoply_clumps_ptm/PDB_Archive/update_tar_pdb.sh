# RUN ME from the ~/Git/PDB_Archive/ directory

## update local pdb database
rsync -rlpt -v -z --delete --port=33444 \
rsync.rcsb.org::ftp_data/structures/divided/pdb/ pdbs/

# ## update local pdb from snapshot
# rsync -avy snapshots.pdbj.org::20250101/pub/pdb/data/structures/divided/pdb/ pdbs_2025/

# ensure the .DS_Store file is deleted
rm pdbs/.DS_Store




### Create Tar-Files:

## tar a folder for every letter
pdb_dir=pdbs
# pdb_dir=pdbs_2025

pdb_tars_dir=$pdb_dir'_tars'
mkdir $pdb_tars_dir
for prefix in $(ls $pdb_dir | awk -F'/' '{print substr($1,1,1)}' | uniq)
do
	tar -cf - --strip-components=1 $pdb_dir/$prefix* | pv -s $(($(du -skc $pdb_dir/$prefix* | tail -n 1 | awk '{print $1}') * 1024)) > tmp_$prefix.tar
	mv tmp_$prefix.tar $pdb_tars_dir/pdbs_$prefix.tar # overwrite old file once move is complete
done
gsutil -m cp -r $pdb_tars_dir/ gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/ # sync to google bucket



# create a testing directory, with only a few full PDB folders
pdb_tars_dir=$pdb_dir'_tars_test'
mkdir $pdb_tars_dir
test_prefix=(a b c d e)
for prefix in $(ls $pdb_dir | awk -F'/' '{print substr($1,1,1)}' | uniq)
do
	echo $prefix
	if printf "%s\n" "${test_prefix[@]}" | grep -q "^$prefix$"; then
		tar -cf - --strip-components=1 $pdb_dir/$prefix* | pv -s $(($(du -skc $pdb_dir/$prefix* | tail -n 1 | awk '{print $1}') * 1024)) > tmp_$prefix.tar
		mv tmp_$prefix.tar $pdb_tars_dir/pdbs_$prefix.tar # overwrite old file once move is complete
	else
		echo "Skipping this directory"
		tar -cf tmp_$prefix.tar -T /dev/null
		mv tmp_$prefix.tar $pdb_tars_dir/pdbs_$prefix.tar # overwrite old file once move is complete
	fi
done
gsutil -m cp -r $pdb_tars_dir/ gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/ # sync to google bucket


## copy over reference files

wget -P . https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip -c uniprot_sprot.fasta.gz > uniprot_sprot.fasta
rm uniprot_sprot.fasta.gz
gsutil -m cp uniprot_sprot.fasta gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/reference_files/

# NOTE: frozen version is encrypted
# wget -P . https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2025_01/knowledgebase/uniprot_sprot-only2025_01.tar.gz
# gunzip -c uniprot_sprot-only2025_01.tar.gz > uniprot_sprot_202501.fasta
# rm uniprot_sprot_202501.fasta.gz
# gsutil -m cp uniprot_sprot_202501.fasta gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/reference_files/



wget -P . ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
gunzip -c pdb_chain_uniprot.tsv.gz > pdb_chain_uniprot.tsv
rm pdb_chain_uniprot.tsv.gz
gsutil -m cp pdb_chain_uniprot.tsv gs://fc-385e9b4e-43ff-44b3-8cf7-036a2a96d102/reference_files/


FROM broadcptacdev/panoply_clumps_ptm:latest
MAINTAINER proteogenomics@broadinstitute.org

# Install BLAST
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.16.0/ncbi-blast-2.16.0+-x64-linux.tar.gz
RUN tar -xf ncbi-blast-2.16.0+-x64-linux.tar.gz; rm ncbi-blast-2.16.0+-x64-linux.tar.gz
ENV PATH="$PATH:/ncbi-blast-2.16.0+/bin"


# Download Reference Files
#RUN wget -P /reference_files/ https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
#RUN gunzip -c /reference_files/uniprot_sprot.fasta.gz > /reference_files/uniprot_sprot.fasta
#RUN rm /reference_files/uniprot_sprot.fasta.gz

#RUN wget -P /reference_files/ ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz
#RUN gunzip -c /reference_files/pdb_chain_uniprot.tsv.gz > /reference_files/pdb_chain_uniprot.tsv
#RUN rm /reference_files/pdb_chain_uniprot.tsv.gz


COPY src /prot/proteomics/Projects/PGDAC/src



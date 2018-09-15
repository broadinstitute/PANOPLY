
# To run these WDL uploads to FireCloud, start with
#  cd ~/proteomics/PGDAC/wdl
#  docker run --rm -it -v "$HOME"/.config:/.config -v "$PWD":/working broadinstitute/firecloud-cli bash

# Tasks
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_aggregate_fasta -t Workflow -y "PGDAC FASTA aggregation" tasks/pgdac_aggregate_fasta/pgdac_aggregate_fasta.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_association -t Workflow -y "PGDAC Association analysis" tasks/pgdac_association/pgdac_association.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_cna_correlation -t Workflow -y "PGDAC CNA correlation analysis" tasks/pgdac_cna_correlation/pgdac_cna_correlation.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_cna_correlation_report -t Workflow -y "PGDAC CNA correlation analysis report" tasks/pgdac_cna_correlation_report/pgdac_cna_correlation_report.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_cna_setup -t Workflow -y "PGDAC CNA correlation analysis setup" tasks/pgdac_cna_setup/pgdac_cna_setup.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_cons_clust -t Workflow -y "PGDAC Consensus k-means clustering" tasks/pgdac_cons_clust/pgdac_cons_clust.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_cpdb -t Workflow -y "PGDAC customProDB pipeline" tasks/pgdac_cpdb/pgdac_cpdb.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_harmonize -t Workflow -y "PGDAC Multi-omic data harmonization" tasks/pgdac_harmonize/pgdac_harmonize.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_mo_nmf -t Workflow -y "PGDAC multi-omics NMF clustering" tasks/pgdac_mo_nmf/pgdac_mo_nmf.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_normalize_ms_data -t Workflow -y "PGDAC Data normalization and filtering" tasks/pgdac_normalize_ms_data/pgdac_normalize_ms_data.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_normalize_ms_data_report -t Workflow -y "PGDAC Data normalization and filtering report" tasks/pgdac_normalize_ms_data_report/pgdac_normalize_ms_data_report.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_parse_sm_table -t Workflow -y "PGDAC input parser for Spectrum Mill" tasks/pgdac_parse_sm_table/pgdac_parse_sm_table.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_rna_protein_correlation -t Workflow -y "PGDAC RNA-protein correlation" tasks/pgdac_rna_protein_correlation/pgdac_rna_protein_correlation.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_rna_protein_correlation_report -t Workflow -y "PGDAC RNA-protein correlation report" tasks/pgdac_rna_protein_correlation_report/pgdac_rna_protein_correlation_report.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_sampleqc -t Workflow -y "PGDAC Sample QC" tasks/pgdac_sampleqc/pgdac_sampleqc.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_sampleqc_report -t Workflow -y "PGDAC Sample QC report" tasks/pgdac_sampleqc_report/pgdac_sampleqc_report.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_ssgsea -t Workflow -y "PGDAC single-sample GSEA with parallel execution" tasks/pgdac_ssgsea/pgdac_ssgsea.wdl

# Workflows
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_main -t Workflow -y "PGDAC Main pipeline" workflows/pgdac_main/pgdac_main.wdl
firecloud -u https://api.firecloud.org/api -m push -s broadcptac -n pgdac_cmap_analysis -t Workflow -y "PGDAC CMAP analysis" workflows/pgdac_cmap_analysis/pgdac_cmap_analysis.wdl
firecloud -u https://api.firecloud.org/api -m push -s cptac -n cptac_outlier_analysis -t Workflow -y "CPTAC NYU outlier analysis" workflows/outlier_analysis/outlier_analysis.wdl

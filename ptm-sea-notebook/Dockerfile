FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-r:2.1.3
USER root

RUN R -e 'install.packages(c("pacman", "optparse", "gtools", "verification", "doParallel", "foreach", "magrittr", "yaml", "glue", "dplyr"))'
RUN R -e 'BiocManager::install(c("rhdf5", "cmapR", "Biostrings"))'

RUN apt-get update && apt-get install -y dos2unix && apt-get install -y nano
RUN mkdir -p /ptm-sea/src
RUN mkdir -p /ptm-sea/db/ptmsigdb
RUN dos2unix /ptm-sea/src

RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/ssGSEA2.0.R && dos2unix ssGSEA2.0.R && mv ssGSEA2.0.R /ptm-sea/src
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/parse_yaml_ssgsea.R && dos2unix parse_yaml_ssgsea.R && mv parse_yaml_ssgsea.R /ptm-sea/src
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/io.R && dos2unix io.R && mv io.R /ptm-sea/src
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/src/utils.R && dos2unix utils.R && mv utils.R /ptm-sea/src
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/ssgsea-cli.R && dos2unix ssgsea-cli.R && mv ssgsea-cli.R /ptm-sea/

RUN wget https://raw.githubusercontent.com/broadinstitute/PANOPLY/dev/src/panoply_ssgsea/panoply_ptmsea_functions.R && dos2unix panoply_ptmsea_functions.R && mv panoply_ptmsea_functions.R /ptm-sea/src
RUN wget https://raw.githubusercontent.com/broadinstitute/PANOPLY/dev/src/panoply_ssgsea/parse_yaml_preprocess_gct.R && dos2unix parse_yaml_preprocess_gct.R && mv parse_yaml_preprocess_gct.R /ptm-sea/src
RUN wget https://raw.githubusercontent.com/broadinstitute/PANOPLY/dev/src/panoply_ssgsea/preproc-ssgsea.yaml && dos2unix preproc-ssgsea.yaml && mv preproc-ssgsea.yaml /ptm-sea/src
RUN wget https://raw.githubusercontent.com/broadinstitute/PANOPLY/dev/src/panoply_ssgsea/preprocessGCT.R && dos2unix preprocessGCT.R && mv preprocessGCT.R /ptm-sea/src

# download PTM signature databases
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/ptm.sig.db.all.flanking.human.v1.9.0.gmt && mv ptm.sig.db.all.flanking.human.v1.9.0.gmt /ptm-sea/db/ptmsigdb
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/ptm.sig.db.all.flanking.mouse.v1.9.0.gmt && mv ptm.sig.db.all.flanking.mouse.v1.9.0.gmt /ptm-sea/db/ptmsigdb
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/ptm.sig.db.all.flanking.rat.v1.9.0.gmt && mv ptm.sig.db.all.flanking.rat.v1.9.0.gmt /ptm-sea/db/ptmsigdb
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/ptm.sig.db.all.sitegrpid.human.v1.9.0.gmt && mv ptm.sig.db.all.sitegrpid.human.v1.9.0.gmt /ptm-sea/db/ptmsigdb
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/ptm.sig.db.all.sitegrpid.mouse.v1.9.0.gmt && mv ptm.sig.db.all.sitegrpid.mouse.v1.9.0.gmt /ptm-sea/db/ptmsigdb
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/ptm.sig.db.all.sitegrpid.rat.v1.9.0.gmt && mv ptm.sig.db.all.sitegrpid.rat.v1.9.0.gmt /ptm-sea/db/ptmsigdb
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/ptm.sig.db.all.uniprot.human.v1.9.0.gmt && mv ptm.sig.db.all.uniprot.human.v1.9.0.gmt /ptm-sea/db/ptmsigdb
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/ptm.sig.db.all.uniprot.mouse.v1.9.0.gmt && mv ptm.sig.db.all.uniprot.mouse.v1.9.0.gmt /ptm-sea/db/ptmsigdb
RUN wget https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/ptm.sig.db.all.uniprot.rat.v1.9.0.gmt && mv ptm.sig.db.all.uniprot.rat.v1.9.0.gmt /ptm-sea/db/ptmsigdb

RUN chmod +x /ptm-sea/src/preprocessGCT.R
RUN chmod +x /ptm-sea/ssgsea-cli.R

ADD terra-functions.R /ptm-sea/src/terra-functions.R

USER $USER

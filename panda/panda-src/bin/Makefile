#!/bin/bash

SHELL=bash
TEDMINT=$(PANDA)/bin
LOGm=logfile_mn.txt
LOGb=logfile_bg.txt
LOGm_EXE=>> $(LOGm) 2>&1
LOGb_EXE=>> $(LOGb) 2>&1

GCT_FILES=$(wildcard pipeline-input/*.gct)
GCT_FLOGS=$(patsubst pipeline-input/%.gct, split-log-%.txt, $(GCT_FILES))
CSV_FILES=$(wildcard pipeline-input/*.csv)
CSV_FLOGS=$(patsubst pipeline-input/%.csv, split-log-%.txt, $(CSV_FILES))

#?
#? =====================================================

#? 
#? ** panda-init
#?     . initialises panda workspace 
#?     . creates groups and Terra workspace
.PHONY: panda-init
panda-init:
	@make restart
	@make init
	@make groups
	@make terrainit

#? 
#? ** panda-samples
#?     . splits input files into per-sample files
#?     . uploads them on Google Bucket
#?     . links them to the Terra Data Model 
.PHONY: panda-samples
panda-samples:
	@make sp-split-gct
	@make sp-split-csv
	@make sp-link

#?
#? ** panda-sets
#?     . creates sample sets
#?     . updates Terra Data model
#?     . links it to Google Bucket files
.PHONY: panda-sets
panda-sets:
	@make ss-create
	@make ss-model
	@make ss-aggregates

#? 
#? ** panda-run
#?     . runs the entire panda pipeline
.PHONY: panda-run
panda-run:
	@make panda-init
	@make panda-samples
	@make ss-create
	@make ss-model
	@make ss-aggregates

#?
#? =====================================================

#?
#? ** restart
#?     . keeps only config.yaml and this Makefile
#?     . deletes rest
.PHONY: restart
restart:
	@echo -e "* Restarting process.."
	@mkdir -p safe
	@cp config.yaml safe/.
	@cp Makefile safe/.
	@find . -maxdepth 1 -type f -delete
	@mv safe/* .
	@rm -r */

#?
#? ** init
#?     . parses config.yaml and prepares workspace
.PHONY: init
init: 
	@echo -e "* Preparing workspace.."
	@cp /dev/null $(LOGm)
	@$(SHELL) $(TEDMINT)/init.sh $(LOGm_EXE)
	@echo -e "* Parsing config.yaml.."
	@$(SHELL) $(TEDMINT)/process-input.sh $(LOGm_EXE)

#? 
#? ** groups
#?     . creates groups file using listed columns 
#?     . from the config.yaml file
.PHONY: groups
groups:
	@echo -e "* Making groups.."
	@$(SHELL) $(TEDMINT)/create-groups.sh $(LOGm_EXE)

#?
#? ** terrainit
#?     . initialises terra workspace
.PHONY: terrainit
terrainit:
	@echo -e "* Preparing Terra workspace.."
	@$(SHELL) $(TEDMINT)/terra-init.sh $(LOGm_EXE)

#?
#? =====================================================

#?
#? ** sp-split-gct
#?     . splits GCT input files into per-sample GCTs
.PHONY: sp-split-gct
sp-split-gct: $(GCT_FLOGS)
split-log-%.txt: pipeline-input/%.gct
	@echo -e "* Splitting .GCT files.."
	@$(SHELL) $(TEDMINT)/split.sh -e "gct" -i $< > $@ 2>&1;
	@echo -e "* Copying them to Terra.."
	@$(SHELL) $(TEDMINT)/terra-copy.sh -e "gct" -i $< >> $@ 2>&1;

#?
#? ** sp-split-csv
#?     . splits CSV input files into per-sample CSVs
.PHONY: sp-split-csv
sp-split-csv: $(CSV_FLOGS)
split-log-%.txt: pipeline-input/%.csv
	@echo -e "* Splitting .CSV files.."
	@$(SHELL) $(TEDMINT)/split.sh -e "csv" -i $< > $@ 2>&1;
	@echo -e "* Copying them to Terra.."
	@$(SHELL) $(TEDMINT)/terra-copy.sh -e "csv" -i $< >> $@ 2>&1;

#?
#? ** sp-link
#?     . uploads samples to the Terra Data Model
#?     . links them to Google Bucket Files
.PHONY: sp-link
sp-link:
	@echo -e "* Adding sample data model to Terra.."
	@$(SHELL) $(TEDMINT)/terra-add-samples.sh $(LOGm_EXE)
	@echo -e "* Setting links up for samples data"
	@$(SHELL) $(TEDMINT)/terra-set-links.sh $(LOGm_EXE)

#?
#? =====================================================

#?
#? ** ss-create
#?     . creates sample set memberships
.PHONY: ss-create
ss-create:
	@echo -e "* Creating sample sets memberships.."
	@$(SHELL) $(TEDMINT)/create-sets.sh $(LOGm_EXE)

#?
#? ** ss-model
#?     . updates Terra Data model with sample sets
.PHONY: ss-model
ss-model:
	@echo -e "* Updating sample set TSVs.."
	@$(SHELL) $(TEDMINT)/terra-add-ss.sh

#?
#? ** ss-aggregates
#?     . creates aggregates 
#?     . uploads them to the Google Bucket
#?     . uploads Other attributes and Additional Params
#?     . links all to the Terra Data model
.PHONY: ss-aggregates
ss-aggregates:
	@echo -e "* Getting GCB addresses for split files.."
	@$(SHELL) $(TEDMINT)/get_split_samples.sh $(LOGm_EXE)
	@echo -e "* Creating aggregates.."
	@$(SHELL) $(TEDMINT)/process-aggregates.sh $(LOGm_EXE)

#?
#? =====================================================

#?
#? ** clean
#?     . cleans logfiles and intermediate files
.PHONY: clean
clean:
	@echo -e "* Cleaning current workspace.."
	@rm *.txt
	@rm config.sh

#?
#? ** help
#?     . prints this help message
.PHONY : help
help : Makefile
	@sed -n 's/^#?//p' $<

#?
#? =====================================================

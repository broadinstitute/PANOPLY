### Workspace for direct import and streamlined processing of raw data files from the Proteomic Data Commons

The **Proteomic Data Commons** ([PDC](https://proteomic.datacommons.cancer.gov/pdc/)) is a data repository within the [NCI Cancer Research Data Commons (CRDC)](https://datacommons.cancer.gov/) and provides access to curated and standardized *proteomic data* along with biospecimen and clinical metadata. The goal of this workspace (notebook+workflow) is to faciliate *automated data download and streamlined analysis* of raw proteomic data imported from the PDC using the [`FragPipe`](https://fragpipe.nesvilab.org/) proteomics pipeline.

**Step-by-step** process for importing and processing raw proteomics data files:
1. Start by navigating to the [PDC](https://proteomic.datacommons.cancer.gov/pdc/). 
2. Browse through the proteomics raw data available there and select files to import and process on Terra. Analyzing TMT or iTRAQ raw data using `FragPipe` requires `mzML` files. On the PDC, `mzML` files are listed in the `Processed Mass Spectra (Open Standard)` data category.
3. Export file manifest for the chosen files using the `PFB` button. This will connect to Terra and prompt for a workspace to put locate the `file` manifest table.
4. Copy the `PDC_Direct_Data_Import` notebook to the workspace with the `file` manifest table and run all the code blocks in sequence. Running the code blocks will:
    1. Download all the files listed in the table to the `fragpipe` directory in the workspace bucket.
    2. Organize files into subdirectories--one subdirectory for each TMT/iTRAQ plex, including all fractions for that plex.
    3. Create annotation files for each TMT/iTRAQ plex to provide sample IDs.
    4. Generate a `FragPipe` manifest file for processing all the data files.
5. Once the notebook has been successfully run, create a `FragPipe` [workflow](https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_workflows.html) using the `FragPipe` graphical user interface. Workflows can also be directly [downloaded](https://github.com/Nesvilab/FragPipe/tree/master/MSFragger-GUI/workflows) and/or edited as needed. Upload the workflow to the workspace bucket.
6. Obtain a relevant (fasta) sequence database and upload to the workspace bucket. Sequence databases can be downloaded using the `FragPipe` graphical user interface, or copied from other workspaces.
7. If not already present, import the `panoply_fragpipe_search` method from the [FireCloud Method Repository](https://portal.firecloud.org/?return=terra#methods).
8. Set up inputs to the `panoply_fragpipe_search` workflow:
    1. Point `database` to the sequence database in the workspace bucket (from Step 6).
    2. `files_folder` is the bucket address of the `fragpipe` directory created to host all the downloaded files.
    3. The location of the manifest file `fragpipe-manifest.fp-manifest` created by the notebook should be filled into `fragpipe_manifest`.
    4. The bucket location of the workflow file uploaded in Step 5 goes in the `fragpipe_workflow` input slot.
9. Run the workflow and download the output `zip` files.

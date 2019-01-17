DATA FILES:

  RAW DATA
    Raw data and metadata files used for preprocessing are in the data/ directory.
    This directory includes the proteomics and genomics data.
    The *.ssv file in the data/ directory is the output from Spectrum Mill.
    Genomics data is formatted as gct files, with samples names matching the ids in
    the experimental design file. RNA and CNA should have gene symbols as id's 
    (Name column in gct v1.2/3). Parsed PTM/proteome (see filter.r) will have gene 
    symbols in Description column (gct v1.2/3)

    Other extracted data columns from the Spectrum Mill outputs are contained in the
    parsed-data/ directory. The information in the files should be evident based
    on file name.


  PREPROCESSED DATA
    (in the normalized-data directory)
    All data is in the form of log2 (sample/ref) ratios. 
    They are formatted as gct files.
    (see http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gct)



    For each group, following files are available:
     *-ratio.gct     				collection of relevant columns from SM output
 		     				no normalization, no filtering
     *-ratio-norm.gct  				data after 2-component normalization (see below)
     *-ratio-norm-nosdfilter-NArm.gct		normalized, filtered for ratio count and missing values
     *-ratio-norm-noNA.gct			normalized, filtered for ratio count and SD.
						Any proteins with (any) missing values are removed.
     *-ratio-norm-NArm.gct			normalized, filtered for ratio count, missing values and SD.


  Annotation (identical cls, i.e., class vectors, for all sets of input datasets):
  (see http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/cls)
     *-<annotation>.cls				Class vector for <annotation>

  HARMONIZED DATA
    Harmonized RNA, CNA and PTM/proteome data is contained in the harmonized-data directory.
    The data are harmonized by selecting common rows (genes) and columns (samples).
    Harmonized data tables are in txt format and will have only a GeneSymbol column in addition 
    to the data. These tables can be used as input to NetGestalt of other gene-centric analysis tool.
    
    


NORMALIZATION
  Using 2-component Gaussian mixture model:
   - Find the mode M using kernel density estimation (Gaussian kernel with Shafer-Jones bandwidth)
   - Fit mixture model with mean for both components constrained to be equal to M
   - Normalize (standardize) samples using mean M and smaller std. dev. from mixture model fit




FILTERING
  Missing value filter: Can be missing in at most na.max (see config.r) samples. 
  Ratio count filter: Must have at least 2 peptide ratios (1 for phosphosites) observed in at least (1-na.max) samples.
  SD Filter: SD of feature (protein or phosphosite) across all samples is > sd.filter.threshold (see config.r).




ANALYSIS (not present in data-freeze datasets):
  cna: Copy number correlation analysis between CNA, RNA and PTM/proteome
  rna: RNA-protein correlation
  clustering: NMF consensus clustering with automated determination of # clusters
  association: Association / marker-selection analysis for class vector (provided or included in input)
  
  
  
  
CODE
  R files containing the code used for preprocessing are also included in the appropriate directories.
  However, no attempt has been made to make the code self-contained. There may be other dependent
  files that are missing. Contact manidr@broadinstitute.org.



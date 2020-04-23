task panoply_rna_protein_correlation {
  File inputData
  File rnaExpr
  String type
  String standalone
  String? subType
  File? params
  String? analysisDir
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_rna_protein_correlation-output.tar"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh RNAcorr \
                  -i ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -rna ${rnaExpr} \
                  -o ${outFile} \
                  ${"-m " + subType} \
                  ${"-p " + params};
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh RNAcorr \
                  -f ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -d ${dataDir} \
                  -rna ${rnaExpr} \
                  -r ${analysisDir} \
                  -o ${outFile} \
                  ${"-m " + subType} \
                  ${"-p " + params};
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker      : "broadcptac/panoply_rna_protein_correlation:1"
    memory      : select_first ([memory, 12]) + "GB"
    disks       : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu         : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email  : "rkothadi@broadinstitute.org"
  }
}

task panoply_ssgsea {
  File input_ds
  ## parameters to create gene-centric or single-site-centric 
  ## GCT files for ssGSEA / PTM-SEA
  String level
  String? id_type
  String? id_type_out
  String? acc_type
  String? seqwin_col
  String? gene_col
  String? SGT_col
  Boolean? loc
  String? mode
  String? mod_res
  String? mod_type
  Boolean? preprocess_gct

  ## ssGSEA / PTM-SEA parameters below
  File gene_set_database
  String output_prefix
  String? sample_norm_type
  String? correl_type
  String? statistic
  String? output_score_type
  Float? weight
  Int? min_overlap
  Int? nperm
  Boolean? global_fdr

  ## VM parameters
  Int? memory
  Int? disk_space
  Int? num_threads
	
  command {
    set -euo pipefail
    ## output prefix
    ##out_pref=${default="out" output_prefix}


    # prepare GCT file
    /home/pgdac/src/preprocessGCT.R -i ${input_ds} \
      -l ${default="ssc" level} \
      -t ${default="sm" id_type} \
      -o ${default="seqwin" id_type_out} \
      -a ${default="refseq" acc_type} \
      --seqwin_column ${default="VMsiteFlanks" seqwin_col} \
      --gene_symbol ${default="geneSymbol" gene_col} \
      -v ${default="subgroupNum" SGT_col} \
      -d ${default=true loc} \
      -m ${default="median" mode} \
      -r "${default="S|T|Y" mod_res}" \
      -p ${default="p" mod_type} \
      -u ${default=true preprocess_gct} \
      -z /home/pgdac/src
	
    # update path to input_ds
    input_ds2=`cat fn.out`
		
    # default value for 'min.overlap'
    #min_overlap2=10
    #if [ ${default="ssc" level} = "ssc" ]; then
    #min_overlap2=5
    #fi

    # run ssgsea/ptm-sea       
    /home/pgdac/ssgsea-cli.R -i $input_ds2 \
      -d ${gene_set_database} \
      -o ${output_prefix} \
      -n ${default="rank" sample_norm_type} \
      -w ${default='0.75' weight} \
      -c ${default="z.score" correl_type} \
      -t ${default="area.under.RES" statistic} \
      -s ${default="NES" output_score_type} \
      -p ${default=1000 nperm} \
      -m ${default=5 min_overlap} \
      --globalfdr ${default=true global_fdr}
		
    # archive result files
    find * -regextype posix-extended \
      -regex "^signature_gct/.*.gct$|^${output_prefix}.*.gct$|^.*.log.txt$|^.*parameters.txt$" \
      -print0 | tar -czvf ${output_prefix}.tar.gz --null -T -
  }

  output {
    # Outputs defined here
    File results="${output_prefix}.tar.gz"
  }

  runtime {
    docker : "broadcptac/panoply_ssgsea:5"
    memory : select_first([memory, 4]) + " GB"
    disks  : "local-disk " + select_first([disk_space, 15]) + " HDD"
    cpu    : select_first([num_threads, 2])
  }

  meta {
    author : "Karsten Krug"
    email  : "karsten@broadinstitute.org"
  }

}

task panoply_rna_protein_correlation_report {
  File tarball
  String label
  String type
  String tmpDir
  Float fdr

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-rna-seq-correlation.r \
      ${tarball} \
      ${label} \
      ${type} \
      ${fdr} \
      ${tmpDir}
  }

  output {
    File report = "rna-corr_" + label + ".html"
  }

  runtime {
    docker : "broadcptac/panoply_rmd:3"
    memory : select_first ([memory, 8]) + "GB"
    disks  : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu    : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email  : "karsten@broadinstitute.org"
  }
}

task panoply_harmonize {
  File inputData
  File rnaExpr
  File cnaExpr
  String type
  String standalone
  String? analysisDir
  String? subType
  File? params

  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String dataDir = "/prot/proteomics/Projects/PGDAC/data"
  String outFile = "panoply_harmonize-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh harmonize \
                  -i ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -d ${dataDir} \
                  -rna ${rnaExpr} \
                  -cna ${cnaExpr} \
                  -o ${outFile} \
                  ${"-m " + subType} \
                  ${"-p " + params};
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh harmonize \
                  -f ${inputData} \
                  -r ${analysisDir} \
                  -t ${type} \
                  -c ${codeDir} \
                  -d ${dataDir} \
                  -rna ${rnaExpr} \
                  -cna ${cnaExpr} \
                  -o ${outFile} \
                  ${"-m " + subType} \
                  ${"-p " + params};
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker      : "broadcptac/panoply_harmonize:1"
    memory      : select_first ([memory, 12]) + "GB"
    disks       : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu         : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email  : "rkothadi@broadinstitute.org"
  }
}

task panoply_sampleqc {
  File tarball   # output from panoply_harmonize
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_sampleqc-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh sampleQC \
      -i ${tarball} \
      -t ${type} \
      -c ${codeDir} \
      -o ${outFile} \
      ${"-m " + subType} \
      ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker      : "broadcptac/panoply_sampleqc:1"
    memory      : select_first ([memory, 12]) + "GB"
    disks       : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu         : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email  : "manidr@broadinstitute.org"
  }
}

task panoply_sampleqc_report {
  File tarball
  String label
  String type
  String tmpDir

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-sample-qc.r \
      ${tarball} \
      ${label} \
      ${type} \
      ${tmpDir}
  }

  output {
    File report = "sample-qc_" + label + ".html"
  }

  runtime {
    docker : "broadcptac/panoply_rmd:3"
    memory : select_first ([memory, 8]) + "GB"
    disks  : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu    : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email  : "karsten@broadinstitute.org"
  }
}

task panoply_cna_setup {
  File tarball   # output from panoply_harmonize
  File? groupsFile
  String type
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_cna_setup-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CNAsetup \
      -i ${tarball} \
      -t ${type} \
      -c ${codeDir} \
      -o ${outFile} \
      ${"-g " + groupsFile} \
      ${"-m " + subType} \
      ${"-p " + params}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker      : "broadcptac/panoply_cna_setup:1"
    memory      : select_first ([memory, 12]) + "GB"
    disks       : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu         : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email  : "manidr@broadinstitute.org"
  }
}

task panoply_cna_correlation {
  File tarball   # output from panoply_cna_setup
  String type
  String? subType
  File? params
  Float fdr_cna_corr = 0.05
  String outFile = "panoply_cna_correlation-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions


  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh CNAcorr \
      -i ${tarball} \
      -t ${type} \
      -o ${outFile} \
      ${"-m " + subType} \
      ${"-p " + params} \
      -z ${fdr_cna_corr}
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker      : "broadcptac/panoply_cna_setup:1"
    memory      : select_first ([memory, 12]) + "GB"
    disks       : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu         : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani"
    email  : "manidr@broadinstitute.org"
  }
}

task panoply_cna_correlation_report {
  File tarball
  String label
  String type
  String tmpDir
  Float fdr

  Int? memory
  Int? disk_space
  Int? num_threads

  command {
    set -euo pipefail
    Rscript /home/pgdac/src/rmd-cna-analysis.r \
      ${tarball} \
      ${label} \
      ${type} \
      ${fdr} \
      ${tmpDir}
  }

  output {
    File report = "cna-analysis_" + label + ".html"
  }

  runtime {
    docker : "broadcptac/panoply_rmd:3"
    memory : select_first ([memory, 8]) + "GB"
    disks  : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu    : select_first ([num_threads, 1]) + ""
  }

  meta {
    author : "Karsten Krug"
    email  : "karsten@broadinstitute.org"
  }
}

task panoply_association {
  File inputData
  String type
  String standalone
  String? analysisDir
  File? groupsFile
  String? subType
  File? params

  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_association-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    if [[ ${standalone} = false ]]; then
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh assoc \
                  -i ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -o ${outFile} \
                  ${"-g " + groupsFile} \
                  ${"-m " + subType} \
                  ${"-p " + params};
    else
      /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh assoc \
                  -f ${inputData} \
                  -t ${type} \
                  -c ${codeDir} \
                  -r ${analysisDir} \
                  -o ${outFile} \
                  -g ${groupsFile} \
                  ${"-m " + subType} \
                  ${"-p " + params}
    fi
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker      : "broadcptac/panoply_association:1"
    memory      : select_first ([memory, 16]) + "GB"
    disks       : "local-disk " + select_first ([disk_space, 40]) + " SSD"
    cpu         : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email  : "rkothadi@broadinstitute.org"
  }
}

task panoply_accumulate {
  File input_tar
  String? output_tar = "panoply_contrasts.tar"
  String module
  String? analysisDir = "input_tarball"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    /prot/proteomics/Projects/PGDAC/src/accumulate.sh \
      -i ${input_tar} \
      -o ${output_tar} \
      -r ${analysisDir} \
      -m ${module}
  }

  output {
    File outputs = "${output_tar}"
    Array[File] list_gct = glob( "${analysisDir}/${module}/contrasts/*.gct" )
  }

  runtime {
    docker      : "broadcptac/panoply_accumulate:1"
    memory      : select_first ([memory, 16]) + "GB"
    disks       : "local-disk " + select_first ([disk_space, 40]) + " SSD"
    cpu         : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email  : "rkothadi@broadinstitute.org"
  }
}

task panoply_cons_clust {
  File tarball   # output from panoply_harmonize or panoply_normalize_ms_data
  String type
  File? groupsFile
  String? subType
  File? params
  String codeDir = "/prot/proteomics/Projects/PGDAC/src"
  String outFile = "panoply_cluster-output.tar"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail
    echo ${type}
    /prot/proteomics/Projects/PGDAC/src/run-pipeline.sh cluster \
      -i ${tarball} \
      -t ${type} \
      -c ${codeDir} \
      -o ${outFile} \
      ${"-m " + subType} \
      ${"-p " + params} \
      ${"-g " + groupsFile}
   
  }

  output {
    File outputs = "${outFile}"
  }

  runtime {
    docker      : "broadcptac/panoply_cons_clust:3"
    memory      : select_first ([memory, 16]) + "GB"
    disks       : "local-disk " + select_first ([disk_space, 40]) + " SSD"
    cpu         : select_first ([num_threads, 8]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "D. R. Mani, Karsten Krug"
    email  : "karsten@broadinstitute.org"
  }
}


task panoply_download
{
  File cons_clust_tar
  File ssgsea_ome_tar
  File ssgsea_rna_tar
  File? ptmsea
  String analysisDir
  String summary_tar = "panoply_main_summary.tar"
  String full_tar = "panoply_main_full.tar"
  Array[File] ssgsea_assoc_tars
  Array[File] ssgsea_clust_tars
  String ssgsea_assoc_dir = "ssgsea_assoc"
  String ssgsea_clust_dir = "ssgsea_clust"

  Int? memory
  Int? disk_space
  Int? num_threads
  Int? num_preemptions

  command {
    set -euo pipefail

    if [ ! -d ${ssgsea_assoc_dir} ]; then
      mkdir ${ssgsea_assoc_dir}
    fi

    if [ ! -d ${ssgsea_clust_dir} ]; then
      mkdir ${ssgsea_clust_dir}
    fi

    index=0
    for file in ${sep=' ' ssgsea_assoc_tars} ; do
      basefilename=$(basename $file)
      index=$((index+1))
      cp $file ${ssgsea_assoc_dir}/$basefilename-$index.tar;
    done

    index=0
    for file in ${sep=' ' ssgsea_clust_tars} ; do
      basefilename=$(basename $file)
      index=$((index+1))
      cp $file ${ssgsea_clust_dir}/$basefilename-$index.tar;
    done

    /prot/proteomics/Projects/PGDAC/src/download.sh \
        -c ${cons_clust_tar} \
        -o ${ssgsea_ome_tar} \
        -r ${ssgsea_rna_tar} \
        -a ${analysisDir} \
        -s ${ssgsea_assoc_dir} \
        -l ${ssgsea_clust_dir} \
        ${"-p" + ptmsea};
  }

  output {
    File summary = "${summary_tar}"
    File full = "${full_tar}"
  }

  runtime {
    docker : "broadcptac/panoply_download:1"
    memory : select_first ([memory, 12]) + "GB"
    disks : "local-disk " + select_first ([disk_space, 20]) + " SSD"
    cpu : select_first ([num_threads, 1]) + ""
    preemptible : select_first ([num_preemptions, 0])
  }

  meta {
    author : "Ramani Kothadia"
    email : "rkothadi@broadinstitute.org"
  }
}

workflow panoply_main {
  File inputData
  File exptDesign
  File rnaExpr
  File rnaExpr_v3
  File cnaExpr
  String analysisDir
  Float corr_fdr
  String dataType
  String? dataSubType
  File? additionalParameters
  File? cna_groups
  File? association_groups
  File? cluster_enrichment_groups
  String standalone = "false"
  Float fdr_cna_corr
  String run_ptmsea

  call panoply_rna_protein_correlation {
    input:
      inputData=inputData,
      type=dataType,
      rnaExpr=rnaExpr,
      analysisDir=analysisDir,
      standalone="true",
      params=additionalParameters
  }

  call panoply_ssgsea as ssgsea_rna {
    input:
      input_ds = rnaExpr_v3
  }

  call panoply_ssgsea as ssgsea_ome {
    input:
      input_ds = inputData
  }
  
  if ( run_ptmsea == "true" ){
    if ( dataType == "phosphoproteome" ){
      call panoply_ssgsea as ptmsea_ome {
        input:
          input_ds = inputData
      }
    }
  } 

  call panoply_rna_protein_correlation_report {
    input:
      tarball=panoply_rna_protein_correlation.outputs,
      label=analysisDir,
      fdr=corr_fdr,
      type=dataType,
      tmpDir="tmp"
  }

  call panoply_harmonize {
    input:
      inputData=panoply_rna_protein_correlation.outputs,
      rnaExpr=rnaExpr,
      cnaExpr=cnaExpr,
      standalone=standalone,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }

  call panoply_sampleqc {
    input:
      tarball=panoply_harmonize.outputs,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }

  call panoply_sampleqc_report {
    input:
      tarball=panoply_sampleqc.outputs,
      type=dataType,
      label=analysisDir,
      tmpDir="tmp"
  }

  call panoply_cna_setup {
    input:
      tarball=panoply_sampleqc.outputs,
      groupsFile=cna_groups,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters
  }

  call panoply_cna_correlation {
    input:
      tarball=panoply_cna_setup.outputs,
      type=dataType,
      subType=dataSubType,
      fdr_cna_corr=fdr_cna_corr,
      params=additionalParameters
  }

  call panoply_cna_correlation_report {
    input:
      tarball=panoply_cna_correlation.outputs,
      type=dataType,
      label=analysisDir,
      fdr=corr_fdr,
      tmpDir="tmp"
  }
  
  call panoply_association {
    input: 
      inputData=panoply_cna_correlation.outputs, 
      groupsFile=association_groups,
      type=dataType,
      subType=dataSubType,
      params=additionalParameters,
      standalone=standalone
  }

  call panoply_accumulate as accumulate_assoc {
    input:
      input_tar=panoply_association.outputs,
      module="association"
  } 

  Array[File] list_gct_assoc = accumulate_assoc.list_gct
  scatter (f in list_gct_assoc){
    call panoply_ssgsea as ssgsea_assoc {
      input:
        input_ds="${f}"
    }
  }

  call panoply_cons_clust {
    input:
      tarball=panoply_association.outputs,
      type=dataType,
      groupsFile=cluster_enrichment_groups,
      subType=dataSubType,
      params=additionalParameters
  }

  call panoply_accumulate as accumulate_clustering {
    input:
      input_tar=panoply_cons_clust.outputs,
      module="clustering"
  }

  Array[File] list_gct_clustering = accumulate_clustering.list_gct
  scatter (f in list_gct_clustering){
    call panoply_ssgsea as ssgsea_clustering {
      input:
        input_ds="${f}"
    }
  }

  call panoply_download {
    input:
      cons_clust_tar=panoply_cons_clust.outputs,
      ssgsea_ome_tar=ssgsea_ome.results,
      ssgsea_rna_tar=ssgsea_rna.results,
      analysisDir=analysisDir,
      ssgsea_assoc_tars=ssgsea_assoc.results,
      ssgsea_clust_tars=ssgsea_clustering.results,
      ptmsea=ptmsea_ome.results
  }

  output {
    File summary_and_ssgsea=panoply_download.summary
    File panoply_full=panoply_download.full
    File rna_corr_report=panoply_rna_protein_correlation_report.report
    File cna_corr_report=panoply_cna_correlation_report.report
    File sample_qc_report=panoply_sampleqc_report.report
  }}

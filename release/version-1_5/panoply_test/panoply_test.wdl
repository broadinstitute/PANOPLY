#
# Copyright (c) 2023 The Broad Institute, Inc. All rights reserved.
#

task panoply_test {
	File? prote_ome
	File? phospho_ome
	File? acetyl_ome
	File? ubiquityl_ome

	Int? memory
	Int? disk_space
	Int? num_threads
	Int? num_preemptions

  Array[Pair[String, File]] ome_pairs =
    [ if defined(prote_ome) then ("proteome", prote_ome) else (),
      if defined(phospho_ome) then ("phosphoproteome", phospho_ome) else (),
      if defined(acetyl_ome) then ("acetylome", acetyl_ome) else (),
      if defined(ubiquityl_ome) then ("ubiquitylome", ubiquityl_ome) else () ]

	
	command {
		echo hello world
	}

	output {
		Array[Pair[String,File]] ome_pairs_full = ome_pairs
	}

	runtime {
		docker : "broadcptac/panoply_common:1_5"
		memory : select_first ([memory, 4]) + "GB"
		disks : "local-disk " + select_first ([disk_space, 20]) + " HDD"
		cpu : select_first ([num_threads, 4]) + ""
		preemptible : select_first ([num_preemptions, 0])
	}

	meta {
		author : "C. M. Williams"
		email : "wcorinne@broadinstitute.org"
	}

}

################################################
## workflow
workflow panoply_test_workflow {
    call panoply_test
}

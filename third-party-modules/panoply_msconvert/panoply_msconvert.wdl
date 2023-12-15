version development
workflow convert_dia { 
  input {
    Directory files_folder
    String output_folder

    Boolean staggered_window = false
    String massError = "10.0"

    File? file_of_files
    Boolean try_one_file=false
    
    Int? num_preemptions
    Int? num_cpus
    Int? ram_gb
    Int? local_disk_gb
  }
  
  String filterparams=""
  
if(staggered_window){
    String filterparams="--filter demultiplex massError="+massError+"ppm" 
    }

call get_array_of_files {
    input:
        files_folder=files_folder,
        file_of_files=file_of_files,
        try_one_file=try_one_file,

        num_cpus=4,
        ram_gb=16,
        local_disk_gb=local_disk_gb,
        num_preemptions=num_preemptions
    }

scatter (file in get_array_of_files.input_files) {
    call convert_raw_to_mzml {
    input:
        file=file,
        filterparams = filterparams,

        num_cpus=4,
        ram_gb=16,
        local_disk_gb=local_disk_gb,
        num_preemptions=num_preemptions
      }
    }
    Array[File] array_of_files = (convert_raw_to_mzml.converted_files)
    call copyfiles{
        input:
        array_of_files = array_of_files,
        output_folder = output_folder,

        num_cpus=4,
        ram_gb=16,
        local_disk_gb=local_disk_gb,
        num_preemptions=num_preemptions
    }
}



task convert_raw_to_mzml {
    input {
        File file
        String? filterparams 

        Int num_cpus=4
        Int ram_gb=16
        Int local_disk_gb=750
        Int num_preemptions=0
    }
    
    command {
        . /etc/profile
        set -euo pipefail
        
        projdir="raw_files"
        mkdir -p $projdir/data

        cd $projdir
        cp ~{file} data/

        if ls data/*.raw &> /dev/null; then
        wine msconvert --zlib --filter "peakPicking vendor msLevel=1-" --filter "zeroSamples removeExtra 1-" ~{filterparams} ~{file} -o data/
        else
        echo "No .raw file is present. Skipping .raw to .mzML conversion."
        fi

        cd ..
        find $projdir/data -name '*.mzML' | xargs mv -t .
    }

    output {
        File converted_files = select_first(glob("*.mzML"))
    }

    runtime {
        docker      : "chambm/pwiz-skyline-i-agree-to-the-vendor-licenses"
        cpuPlatform : "AMD Rome"
        memory      : "${ram_gb}GB"
        bootDiskSizeGb: 128
        disks       : "local-disk ${local_disk_gb} HDD"
        preemptible : num_preemptions
        cpu         : num_cpus
    }
    meta {
        author: "Khoi Pham Munchic"
        email : "kpham@broadinstitute.org"
    }
}


    task get_array_of_files {
    input {
        Directory files_folder
        File? file_of_files
        
        Boolean try_one_file=false
        Int num_cpus=32
        Int ram_gb=128
        Int local_disk_gb=500
        Int num_preemptions=0
    }
    Array[File] files = if defined(file_of_files) then read_lines(select_first([file_of_files])) else []
    
    command {
        . /etc/profile
        set -euo pipefail
        
        # localize files
        if [ -z ${file_of_files} ]
        then
        cp -s -a ~{files_folder} .
        mv $(basename ~{files_folder}) data
        else
        mkdir data
        cp -s ~{sep(' ', files)} data
        fi

        if [ "${try_one_file}" = true ]
        then
        file_to_keep=$(ls -1 "data" | head -1)
        cd data && ls | grep -v $file_to_keep | xargs rm -r && cd ..
        fi

        ls -R
    }

    output {
        Array[File] input_files=glob("data/*.raw")
    }

    runtime {
        docker      : "broadcptacdev/panoply_fragpipe:latest"
        cpuPlatform : "AMD Rome"
        memory      : "${ram_gb}GB"
        disks       : "local-disk ${local_disk_gb} HDD"
        preemptible : num_preemptions
        cpu         : num_cpus
    }

    meta {
        author: "Khoi Pham Munchic"
        email : "kpham@broadinstitute.org"
    }
}


task copyfiles{
    input {
        Array[File] array_of_files
        String output_folder 

        Int num_cpus=32
        Int ram_gb=128
        Int local_disk_gb=500
        Int num_preemptions=0
    }
    command{
        . /etc/profile
        set -euo pipefail
        
        gsutil cp ~{sep(' ', array_of_files)} ~{output_folder}

    }

    runtime{
        docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
        cpuPlatform : "AMD Rome"
        memory: "${ram_gb}GB"
        bootDiskSizeGb: 512
        disks : "local-disk ${local_disk_gb} HDD"
        preemptible : num_preemptions
        cpu: num_cpus
    }
    meta {
        author : "Simone Gohsman"
        email : "gohsmans@broadinstitute.org"
    }
}
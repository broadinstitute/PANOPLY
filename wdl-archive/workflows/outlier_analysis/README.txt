
To use docker image and run outlier analysis:

# pull docker
docker pull huazhou2/outlier:v3

# run with parameters
docker run --name outlier_analysis --rm -v ${path_with_input.tar}:/input -v ${path_for_output}:/output \
    huazhou2/outlier:v3 /prog/outlier.sh /input/<input-tar-file> /output/<output-tar-file> 10


The input tar file should contain a sample_input_<prefix>.txt that has the following info:
 1st column = file name of data file  [[filename]]
 2nd column = datatype  [[string]]
 3rd column = data description  [[string]]
 4th column = color to use when plotting outlier status  [[#RRGGBB]]
 5th column = cutoff to be used for outlier definition  [[numeric]]
 6th column = does original data need to be normalized?   [[True/False]]
 All following columns specify the grouping that you want to be plotted against in outlier grid plot.
 Group names are based on the column annotations in the GCT1.3 input data file
One row with above info for EACH input GCT1.3 file (ex., proteome, phosphoproteome, RNA, CNA, ...)

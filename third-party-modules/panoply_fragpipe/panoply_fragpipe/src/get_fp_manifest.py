import sys
import os

abs_path = sys.argv[1]
data_folder = sys.argv[2]
raw_file_type = sys.argv[3]

files = os.listdir(data_folder)
manifest = [os.path.join(abs_path, file) + "\t\t\t" + raw_file_type for file in files]

with open('generated.fp-manifest', 'w') as filehandle:
    for f in manifest:
        filehandle.write('%s\n' % f)
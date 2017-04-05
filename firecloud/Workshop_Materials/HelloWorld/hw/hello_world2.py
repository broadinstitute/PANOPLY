import sys
from shutil import copyfileobj as copy


f = open ("/data/out.txt", 'w')
if len(sys.argv) == 1:
    f.write ("Hello World")
else:
    x = sys.argv[1]
    with open (x, 'r') as infile:
        copy (infile, f)



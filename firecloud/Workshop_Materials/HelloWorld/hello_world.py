import sys

if len(sys.argv) == 1:
    x = "World"
else:
    x = sys.argv[1]

file = f.open("hello-world-output.txt", 'w')
file.write("Hello %s!" % x)
file.close()


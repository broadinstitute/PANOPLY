#!/bin/bash
task="bogus-task"

# ./setup.sh -t $task # ok
# ./setup.sh -f -t $task # ok
# ./setup.sh -f # ok
# ./setup.sh -f -t $task -p # ok
# ./setup.sh -t $task -s # ok
# ./setup.sh -t $task -s -f # ok
# ./setup.sh -s # ok
# ./setup.sh -t $task -m $task:latest # ok
# ./setup.sh -t $task -b # ok
# ./setup.sh -t $task -y # ok
# ./setup.sh -t $task -n broadcptac -b # ok
# ./setup.sh -t $task -b -n broadcptac # ok
# ./setup.sh -t $task -u # ok
# ./setup.sh -t $task -u -n broadcptac # ok
# ./setup.sh -t $task -b -u -n broadcptac -g 2 # ok
# ./setup.sh -z # ok

# ./update.sh -t pgdac_common -n broadcptac -g 2 # ok
# ./update.sh -t pgdac_utils -n broadcptac -g 3 # ok

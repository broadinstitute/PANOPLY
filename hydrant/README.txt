
To build dockers for tasks, here is a step-by-step:
0. This operation requires docker (and hence cannot be run on Broad UGER)
   Ensure that hydrant is installed. If not, execute the following:
   % pip install --upgrade git+https://github.com/broadinstitute/HydrantFC.git
1. Login to docker hub so that access to broadcptac is enabled:
   % docker login -u <docker-hub-user>
2. Execute commands from the PGDAC/hydrant directory
3. Invoke setup_workflows_primary.sh (use -h for help)
   Following options are usually required:
   -s <sub_dir>     <sub_dir> = tasks
   -n <wf_name>     <wf_name> is subdirectory under PGDAC/src
                    specify name only (no path needed)
   -d <docker_ns>   <docker_ns> = broadcptac
                    <docker_hub_user> (step 1) much have write access
   -t <docker_tag>  is needed; else docker build/push fails
   -b               to build and push docker


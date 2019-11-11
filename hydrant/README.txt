
To add a new module to the pipeline:
0. This operation requires docker (and hence cannot be run on Broad UGER)
   Ensure that hydrant is installed. If not, execute the following:
   % pip install --upgrade git+https://github.com/broadinstitute/HydrantFC.git
1. Login to docker hub so that access to broadcptac is enabled:
   % docker login -u <docker-hub-user>
2. Create a folder src/pgdac_new_module
3. Create or copy source files to src/pgdac_new_module
4. Go to hydrant/
5. Using setup.sh, create a hydrant workspace for the new module 
   (-h has the documentation on how to use the script)
   ./setup.sh -t pgdac_new_module -f -s -y -n broadcptac -g 1 -b
   -- for creating a new workspace, and copying source files to the 
      docker directory for build operation, build and push docker.
   -- Edit .wdl manually
6. Edit the pgdac_new_module.wdl manually or point setup.sh to the WDL 
   you have already created. 
   ./setup.sh -t pgdac_new_module -f -s -w /path/to/wdl/ -n broadcptac -g 1 -b
   -- for creating a new workspace, copying source files to the docker directory, 
      copying a ready wdl file to the hydrant workspace, build and push docker to 
      docker hub.
Notes:
Build a separate docker called pgdac_new_module (using setup.sh) which will be
built on top of pgdac_common docker (as of now). It thus will technically have
run-pipeline.sh in the src directory but won't need it. (Note:- This is because
pgdac_main's original modules need run-pipeline.sh and the new ones do not.
Post run-pipeline.sh disintegration, we will be able to bypass the pgdac_common
docker altogether making it obsolete, removing run-pipeline.sh from all
dockers. Even so, update.sh takes care of all the dependencies. On update of a
parent docker it determines which children dockers need to be updated, and then
updates them all.) 
All the steps with respect to hydrant can be carried out using setup.sh and update.sh.

For updating exisiting modules, procedure is similar -- use the help option
to provide the appropriate arguments to setup.sh


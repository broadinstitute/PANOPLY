### PANOPLY without Terra: Running PANOPLY on your local machine

PANOPLY can be run without Terra, in your local computer. This would require the ability to run command line programs, edit files and navigate directories from a terminal window. 

In order to run components of PANOPLY on your local machine, follow these steps to for specific tasks or workflows:

1. Clone the appropriate release branch of the PANOPLY repository or download the zip from GitHub. Current release can be found [here](https://github.com/broadinstitute/PANOPLY/tree/release-1_0). In the following `<PANOPLY>` will refer to the root directory of the cloned repository.
2. Install `hydrant` from [GitHub](https://github.com/broadinstitute/HydrantFC), using these [installation instructions](https://github.com/broadinstitute/HydrantFC/wiki/Installation)
3. Navigate to the ``directory.
4. Find the task or workflow to run:
    - To run a `<workflow>` (pipeline), `cd` to the `<PANOPLY>/hydrant/workflows/<workflow>` directory. For example, to run the `panoply_unified_workflow`, execute `cd <PANOPLY>/hydrant/workflows/panoply_unified_workflow`.
    - To run a `<task>` (method/module), `cd` to the `<PANOPLY>/hydrant/tasks/<task>` directory. For example, to run the consensus clustering task `panoply_cons_clust`, execute `cd <PANOPLY>/hydrant/tasks/panoply_cons_clust`.
5. Once in the correct directory (with a `wdl` file), create a subdirectory names `tests`, and validate the `wdl` using hydrant:
```hydrant validate```
This will create a new file names `inputs.json` in the `tests` directory.
6. Edit the `inputs.json` file to fill in required parameter values, and any optional parameters you want to set/change. For any files that are provided as input, enter the *complete (absolute) paths*---relative paths result in unexpected errors. The `yaml` file with all the default parameter settings can be found at `<PANOPLY>/src/panoply_common/master-parameters.yaml`.
7. Once the `inputs.json` file is ready, execute the task/workflow using:
```hydrant test```
8. Once the test is complete, results and output files can be found in the the `cromwell-executions` directory.


PANOPLY tasks and workflow can also be run locally **without using `hydrant`**. Follow these steps:

* Make sure you download and store [cromwell](https://github.com/broadinstitute/cromwell/releases) and [wdltool](https://github.com/broadinstitute/wdltool/releases/tag/0.14) in a folder where you want to store the executions and results of your pipelines. *NOTE: Make sure you pass absolute file paths everywhere and not relative paths to your current working directory.* Failing to do so might result in unexpected errors. 

Following commands replace calls to `hydrant`:

1. `java -jar wdltool.jar validate <WDL File>` you can perform full validation of the WDL file including syntax and semantic checking.
2. `java -jar wdltool.jar inputs <WDL File> > <your_inputs>.json` to print a `.json` skeleton file of the inputs needed for this workflow. Fill in the values in this `.json` document and pass it in to the `run` subcommand.
3. `java -jar cromwell.jar run [options] <your_wdl_file> -i <your_inputs>.json` to run the workflow through the cromwell engine and print out the outputs in `.json` format.

Steps (1) and (2) collective perform the function of `hydrant validate` and step (3) is equivalent to `hydrant test`.



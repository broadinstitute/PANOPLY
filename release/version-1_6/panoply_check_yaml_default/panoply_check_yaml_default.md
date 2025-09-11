# ```panoply_check_yaml_default```

## Description
This internal support module is used to check a parameter value within the provided YAML, returning a boolean TRUE/FALSE value.  The module assumes that the parameter is a boolean, and that the parameter is in the first level of nesting. It is used to determine whether a module should be run or skipped in a workflow setting; for instance, [panoply_main_internal](./Pipelines%3A-panoply_main_internal) uses this module to check the notebook-set parameter `run_ptmsea`, to determine if [PTM-SEA](./Data-Analysis-Modules%3A-panoply_ssgsea) should be run on phosphoproteome data.

## Input

### Required inputs:

* ```yaml```: (`.yaml` file) master-parameters.yaml
* ```param_lookup```:  (String) name of parameter in YAML to be searched.

## Optional inputs:

* ```param```:  (String) Current (Terra set) parameter value. If set, the YAML value will be ignored.


## Output

* ```param_boolean```: Boolean value of YAML parameter (or `param` override)
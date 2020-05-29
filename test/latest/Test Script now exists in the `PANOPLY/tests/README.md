Test Script now exists in the `PANOPLY/tests/latest` folder that compares a gold-standard output tarball with a test tarball. 

### Function
- Checks for missing `png`, `jpeg`,  and `pdf` files in addition to other "data" files
- Checks for missing sub-directories
- Compares all `GCT`s and `CSV`s unless specified with other regular expression rules to narrow down the list of files to compare
- Uses a `tolerance` argument` that can be set while calling the script
- Logs messages for failed comparisons, missing files, and passed files

### Details
- Inputs:
    + `gold.tar`: gold standard tarball
    + `test.tar`: test tarball
    + `tolerance`: threhsold for comparing numerical matrices; by default set to `1.490116e-08` or `.Machine$double.eps ^ 0.5`
    + `file.list`: a YAML file listing patterns
- More on the YAML file:
    + keys should be sub-directories inside the tarball (NOTE: if sub-directory name is `harmonized-data` use key `make.names( "harmonized-data" )` as the key, which in this case would be `harmonized.data`
    + each sub-directory should have a non-nested list of patterns
    + patterns should be regular expressions for file matching with `list.files()`
    + patterns can also be of the form `"!<pattern>"` (ex. `"!sample-info.csv"` or `"!*.tsv"` )
    + by default, this script will compare all CSV files assuming that first row will be its header
    + by default, this script will compare all GCT files' `@mat` attributes   
- Report will be stored as `test-report.txt`


**NOTE** that, this script does one test for all the sub-directories. If you want to test for just one module, this should still work as long as you give it correct inputs. Other directories will simply be listed as missing. 

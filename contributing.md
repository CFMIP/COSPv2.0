# Working practices for COSP development
This document outlines the recommended working practices for COSP development. The sections below describe the workflow that any developer is expected to follow when developing changes to COSP.

## Open an issue
An issue (Github terminology) is an editable webpage that is used to provide a traceable documentation of every change: bugfix, optimization, new capability, etc. You should create an issue when you are planning a change. This issue will describe the purpose of the change. The proposed changes should be discussed with the PMC, and with the code owners of the part of COSP that is going to be modified.
It is the developers responsibility to contact the PMC at the start of this processes (the PMC will not be monitoring the issues) to ensure that proposed work is in keeping with COSP objectives and, that once developed and tested, is likely to be added to the trunk.

## Branch from master
In order to develop a change, you will need your own forked copy of the COSP repository. Once you have done that, you will create a new branch in which you will develop your new changes. 

## Run regression tests
Sometimes, using a different compiler or running in a different architecture in which the reference files have been created can produce differences that will show in the regression tests. In order to catch these situations, pleaser run the regression tests locally before doing any development. If differences are reported, please document these in the issue, and update your reference files. 

## Develop your change
You should make regular commits to your branch to provide a history of your work that others can follow. This will help them to understand what you've done and why.

## Test your change
COSP provides you with a battery of regression tests to check that your changes don’t break the code. Test your changes regularly, not only at the end of the development process. This will help identify problems early.

Currently, we have the follwing regression tests, which must be run from the `driver/run/` directory. They produce outputs that are saved in `driver/data/outputs/UKMO/`.

**Table 1.** Regression tests.

Test # | Command | Input file | Output file | KGO | Description
------------ | ------------ | ------------- | ------------- | ------------- | -------------
1 | cosp2test cosp2_input_nl.txt | cosp2_input_um.nc | cosp2_output_um.nc | cosp2_output_um.gfortran.kgo.nc | Basic test (~150 gridpoints).
2 | cosp2test cosp2_input_nl.um_global.txt | cosp2_input.um_global.nc | cosp2_output.um_global.nc | cosp2_output.um_global.gfortran.kgo.nc | Low-res global model snapshot.

The data files for Test 1 are small and live in the github repository, and therefore distributed with the the code. Other data files are stored in Google Drive and can be downloaded by running `download_test_data.sh` from within the `driver/` directory. The script `download_test_data.sh` checks the downloaded data files against md5 sums in `data/outputs/UKMO/`, to make sure that you are using the correct files.

In `driver/`, the script `compare_to_kgo.py` can be used to test your ouputs against the references files. For instance, for the first test in the table above:

`python compare_to_kgo.py data/outputs/UKMO/cosp2_output_um.gfortran.kgo.nc data/outputs/UKMO/cosp2_output_um.nc`
    
The script accepts thresholds for absolute and relative tolerances, named atol and rtol respectively. By default the script will report all differences, i.e. --atol=0.0 --rtol=0.0.

Any modification that reduces the performance and it is not protected by a logical switch (e.g. new simulator or diagnostic) will have to be strongly motivated.

## Document your change
Add meaningful comments to the source code, and modify the documentation if necessary.

## Code review
Once you are happy with the changes and documentation, and all the necessary tests have passed, you can initiate a _pull request_. This will trigger a battery of regression tests that are run in the github servers. These tests are the ones listed in Table 1, but run with a driver program built with different compilers. The outputs of these tests are recorded in the Actions sections of the gthub repository. If some of your changes do not pass the necessary tests, then go back and figure out why and fix the problem so that they pass the tests.

The developer should assess the impact of the change and the need for a review and get agreement from the PMC. The review process is open to community members outside the PMC. If a review is considered necessary, the developer should find someone to do the review, not necessarily a PMC member.

Typically, the reviewer will suggest improvements and modifications that will lead to iterations in the modified code.  Once the reviewer is satisfied, they will approve the changes and the branch will be ready for merging.

**Do you changes change the results?**

If answers change, provide a detailed explanation why the code changes modify the results. This will trigger discussion with the PMC about your changes.
This process is facilitated by the Continous Integration (CI) tests. When the comparison against the KGO fails, the script `plot_test_outputs.py` is run and plots for the new output file are produced and uploaded to the _Artifact_ created by the CI action. Currently, `plot_test_outputs.py` can only produce plots for 2D output files (like Test 2). If the new results are accepted, you will need to change the md5 sums of the relevant files, and a member of the PMC will upload the new files to Google Drive.

## Merge to the master branch (trunk)
Now that your changes have been reviewed, tested, and approved, a member of the PMC will merge your code into the master branch.

## Close the issue
Once the code has been merged, the issue that documents the changes can be closed.

## Long term commitment
For major changes (e.g. new simulator), the developer will become the owner of that section. Section owners are expected to engage in a long-term commitment to help maintaining the code.

# cbn_analysis - an automated command line tool for visualizing and measuring ellipticity of cucurbituril host/guest structures

## How to install
**This package is not yet pip installable, but once it is available, the command will be a simple**

_pip install cbn_analysis_

**Until then, the installation of the package can be done by the following:**
* Clone this git repository to your local machine
* In the cloned repository, install the cbn_analysis package through _python -m pip install . vv_ 
* Next, install the required dependencies by running _pip install -r requirements.txt_

## How to run the analysis
### The cbn_analysis package takes a single command line argument: an xyz file containing atom coordinates, or a folder of xyz files containing atom coordinates

**To run the analysis package on a single file, navigate to the directory with the cbn_analysis script and execute it:**
(NOTE: if there are spaces in your file name, you will need quotes around the file name)
![](images/single_file)


**The same convention can be used to execute the package on a folder of xyz files**
![](images/folder_test)



## How to interpret the output

[Interactive version found here](https://plotly.com/~Mshavlik/63/)

![](images/testing_cbn_interactive.png)


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
![](images/single_file.png)


**The same convention can be used to execute the package on a folder of xyz files**
![](images/folder_test.png)



## How to interpret the output

### Each xyz file analyzed with this package will produce at least 3 pieces of data: 
* 1 interactive 3D scatter plot with all carbons in the CBn structures visualized 
* 1 interactive 3D scatter plot with only the central carbon ring of the CBn structures visualized
* 1 spreadsheet with all of the carbons, their positions, distance to the centroid of the structure, and the measured ellipticity

**In the case of host CBn structures with an internally situated guest structure, an additional 3D graph will be produced with the guests visualized in the structures (see below)**

![](images/testing_cbn_interactive.png)

To see the interactive version of this plot that gets generated from the script, [Click Here](https://plotly.com/~Mshavlik/63/)

**The user can visually see the ellipticity in the structues and compare them to the measured values in the spreadsheets:**                 
_CB7 structures from above with guests removed_
![](images/CB7_circular.png)

Calculated ellipticity for each structure on a scale of circular (0) to linear (1)\                                                         
![](images/circular_ellipticity.png)


_CB10 wide structures with guests removed_
![](images/ellipsoid_example.png)

Calculated ellipticity for each structure on a scale of circular (0) to linear (1)\                                                         ![](images/ellipse_ellipticity.png)



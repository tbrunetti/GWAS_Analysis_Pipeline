# GWAS_Analysis_Pipeline

## Overview and Purpose
------------------------
An automated pipeline to analyze GWAS data after initial round of QC/filtering has been performed.

<p align="center">
<img src="https://github.com/tbrunetti/GWAS_Analysis_Pipeline/blob/master/how_to_run_flowchat.png" />
</p>


## Software Requirements
------------------------
The following are the minimum software requirements:
* Python version 2.7 (https://www.python.org/)
* R version 3.2 or better (https://cran.r-project.org/)
* PLINK version 1.9 or better (https://www.cog-genomics.org/plink2)
* KING software package (http://people.virginia.edu/~wc9c/KING/manual.html)
* chunkypipes (http://chunky-pipes.readthedocs.io/en/stable/getting_started.html)
* virtualenv -- ONLY required if admin rights are not granted (https://virtualenv.pypa.io/en/stable/) 

__*--Software Requirements that can be installed automatically--*__  
The following list of Python libraries are required but the pipeline can automatically install them if pip is available:
  * SciPy stack, in particular the following packages: (https://scipy.org/)
    * pandas
    * numpy
    * matplotlib
  * Statistics (https://pypi.python.org/pypi/statistics)
  * Pillow (https://pypi.python.org/pypi/Pillow/3.4.2)
  * pyFPDF (http://pyfpdf.readthedocs.io/en/latest/index.html)
  * pyPDF2 (https://pypi.python.org/pypi/PyPDF2/1.26.0)


## User Generated File Requirements
-----------------------------------


## Installing and Running Pipeline on HPC or without sudo privileges
---------------------------------------------------------------------
1. Install, Create, Activate Virutal Environment
2. Clone Repository

### Installation, Creation, and Activation of Virtual Environment (installAtion and creation only needs to ever be done once)
------------------------------------------------------------------------------------------------------------------------------
For users running the pipeline on a HPC system or without administrative privileges, a virtual environment must be created in order to utilize the pipeline to its maximum capabilities and efficiency.  The reason being, certain packages and dependencies can be installed using the pipeline to provide the user with an easy-to-use experience.  Many of these packages will install on your system globally, which cannot be done on shared distributed file systems such as HPC or if administrative rights are not granted on the computer you are using.  To bypass this, you will need to create a Python Virtual Environment.  Below are the steps for installing and creating this environment.  NOTE:  THIS ONLY EVER NEEDS TO BE DONE ONE TIME!

1.  Download and Install virtualenv
------------------------------------
For the full documentation of of virtualenv please refer to the following website:  https://virtualenv.pypa.io/en/stable/installation/  Typing in the following commands below will download virtualenv from the Web and create a virtual environment call myVE under the virtualenv-15.0.0 directory.  The name myVE can be changed to any full path file name of your choosing so that the virtual environment does not need to be listed in the virtualenv-15.0.0 directory.

```
$ curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-15.0.0.tar.gz
$ tar xvfz virtualenv-15.0.0.tar.gz
$ cd virtualenv-15.0.0
$ python virtualenv.py myVE
```
Now you have created a virtual environment call myVE in the directory virtualenv-15.0.0.  Next you need to activate the environment.  In order to activate your virutal environment type in the commmand below:
```
$ source bin/activate
```
Upon successful activation, your commmand prompt should now read something like this:
```
(myVE)user@myaddress $  
```
where the parantheses is the name of the virtual environment you just created and appears to the far left to let you know you have activated and are working in the virtual environment.  

Now anytime you would like to work in the virutal environment make you use
```
$ source bin/activate
```
and see the parenthesis and it is all set up!  Note, however, although you every only need to install the virtual environment once, each time you log back into your system you will need to activate the environment if you want to use it.  Inorder to logout or deactivate the virtual environment just type the following into your command prompt:
```
$ deactivate
```
Now you are working in your regular environment.  All the data you generated in the virtual environment is saved and you can access it like any other directory or folder in your regular working environment.

It is also **critical** on HPC systems to load the python module **PRIOR TO ACTIVATION** of your virtual environment.

### Modifying the slurm scirpt
------------------------------


## Installing and Running Pipeline with on personal system with sudo privileges
-------------------------------------------------------------------------------
1. Clone Repository

##  Output and Deliverables
---------------------------

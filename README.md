# GWAS_Analysis_Pipeline



## Running Pipeline on HPC or without sudo privileges
------------------------------------------------------
1. Install, Create, Acitvate Virutal Environment
2. Clone Repository

### Installation, Creation, and Activation of Virtual Environment (installtion and creation only needs to ever be done once)
-----------------------------------------------------------------------------------------------
For users running the pipeline on a HPC system or without administrative privileges, a virtual environment must be created in order to utilize the pipeline to its maximum capabilities and efficiency.  The reason being, certain packages and dependencies can be installed using the pipeline to provide the user with an easy-to-use experience.  Many of these packages will install on your system globally, which cannot be done on shared distributed file systems such as HPC or if administrative rights are not granted on the computer you are using.  To bypass this, you will need to create a Python Virtual Environment.  Below are the steps for installing and creating this environment.  NOTE:  THIS ONLY EVER NEEDS TO BE DONE ONE TIME!

1.  Download and Install virtualenv
------------------------------------
For the full documentation of of virtualenv please refer to the following website:  https://virtualenv.pypa.io/en/stable/installation/  By typing in the following commands below, this will download virtualenv from the Web and create a virtual environment call myVE under the virtualenv-X.X directory.  The name myVE can be changed to any full path file name of your choosing.

```
$ curl -O https://pypi.python.org/packages/source/v/virtualenv/virtualenv-X.X.tar.gz
$ tar xvfz virtualenv-X.X.tar.gz
$ cd virtualenv-X.X
$ python virtualenv.py myVE
```
Now you have created a virtual environment call myVE in the directory virtualenv-X.X.  Next you need to activate the environment.  In order to activate your virutal environment type in the commmand below:
```
source bin/activate
```
Upon successful activation, your commmand prompt should now read something like this:
```
(myVE)user@myaddress $  
```
where the parantheses with the name of the virtual environment you just created appears to the far left to let you know you have activated and are working in the virtual environment.  

Now anytime you would like to work in the virutal environment make you use
```
source bin/activate
```
and see the parenthesis and it is all set up!  Note, however, although you every only need to install the virtual environment once, each time you log back into your system you will need to activate the environment if you want to use it.  Inorder to logout or deactivate the virtual environment just type the following into your command prompt:
```
deactivate
```
Now you are working in your regular environment.  All the data you generated in the virtual environment is saved and you can access it like any other directory or folder in your regular working environment.


### Modifying the slurm scirpt
------------------------------


## Running Pipeline with on personal system with sudo privileges
-----------------------------------------------------------------
1. Clone Repository

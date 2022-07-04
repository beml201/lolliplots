# lolliplots
A python tool to create lolliplot visualisations

## Setup
You will need an anaconda installation to setup the environemnt.yml file and use the script.
This may be available through a module system, where you should load it into your instance before starting.
eg:
``` module load Anaconda ```

Once *anaconda* is activated, the environment should be setup using the supplied yaml file (*liftover_environment.yml*) using the following code:

``` conda env create -f lolliplots_environment.yml ```

This only needs to be done once.

You can check that the environment is installed correctly using:

``` conda env list ```

## After Setup
Ensure Anaconda is loaded:

``` module load Anaconda ```

Activate the environment:

``` conda activate lolliplots ```

The python file can then be used. For help on the command line inputs use (within the directory the python file is in):

NOTE: python version may need to be specified
``` python lolliplots_example.py ```

Deactivate the environment with
```conda deactivate ```
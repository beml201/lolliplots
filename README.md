# Lolliplots
A python tool to create lollipop plots (lolliplots) for coding region association results

## Setup
You will need an anaconda installation to setup the *lolliplots_environment.yml* file and use the script.
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
***or***
``` source activate lolliplots ```


The example python file can then be used. For help on the command line inputs use (within the directory the python file is in):

NOTE: python version may need to be specified
``` python lolliplots_example.py ```

Deactivate the environment with
```conda deactivate ```

## Loading it into your own script
Make sure to activate the lolliplots environment:
```conda activate lolliplots ```
Lolliplots is not currently a conda-available package, and so needs to be imported by adding the *lolliplots.py* file to your current system. This can be done with the following:
```
import sys
sys.path.append("lolliplot_folder")
import lolliplots as lp
```

Replacing **lolliplot_folder** with the folder you have downloaded your *lolliplots.py* script in

You should then be able to use the functions as normal.

Exon location data will need to be downloaded only once. This can be run using ```lp.get_exon_data() ```

This will open a prompt to download the biomart data.
```lp.read_exon_locs()``` should be given the file that was downloaded using ```lp.get_exon_data()```.

This data should then be filtered to the desired transcript before plotting:
``` selected_exons_df = exons.loc[exons['Transcript stable ID'] == "transcript_id"] ```
Where **transcript_id** is the name of your desired transcript

Use ``` new_assoc, new_exon = lp.reduce_gaps(assoc, selected_exons_df)``` to remove introns from the association data and exon location data before passing the new data to ``` lp.lolliplot_raw() ```

## Functions

### lp.get_exon_data():
Arguments: **None**

Run this command (No need to put into a variable) to download the exon location data.
Downloads the data to a new folder in the current working directory: *ensembl_exon_positions/b38_downloaded*
This only need to be run once. In subsequent runs, this should be loaded in using ```lp.read_exon_locs()``` through your own hard coded location.
User inputs will be requested.

### lp.read_exon_locs(file):
|Arguments| About|
|---|---|
|**file** | string: the file location of the downloaded exon data (eg "ensembl_exon_positions/b38_downloaded") |

Function reads in the downloaded file data using pandas and assigns the appropriate headers to work with the default options of ``` lp.lolliplot_raw() ```
Column names are as follows:
- Gene stable ID
- Transcript stable ID
- Exon start (bp)
- Exon stop (bp)
- Gene start (bp)
- Gene end (bp)
- Gene name
- Exon Name


### lp.reduce_gaps(gene_results, exon_info)
|Arguments|About|Default|
|---|---|---|
|**gene_results** | pd.DataFrame: A pandas dataframe of the association results to plot| - |
|**exon_info** | pd.DataFrame: A pandas dataframe of the filtered exon locations| - |
|new_gap| numeric: New gap size between the exons | 10|
|ex_start_col | string: Name of column of exon_info of exome start locations | 'Exon start (bp)'|
|ex_stop_col | string: Name of column of exon_info of exome stop locations | 'Exon stop (bp)'|
|gene_pos_col | string: Name of column of gene_results of the variant positions | 'GENPOS'|

Reduce the spaces between the coding regions to a smaller gap (designated by **new_gap**). Returns two items, the modified **gene_results** and the modified **exon_info**.
The resulting dataframes can be used to make the lolliplots with the intronic regions removed.

> Note: the exon_info should be filtered to your desired transcript before inputting into ```lp.reduce_gaps()``` or ```lp.lolliplot_raw()```. This can be done in pandas using the ```.loc``` attribute eg: ```df = df.loc[df['Transcript ID'] == TRANSCRIPT_ID, :]```

### lp.lolliplot_raw(results, exon_info, title)
|Arguments|About|Default|
|---|---|---|
|**results**| pd.DataFrame: A pandas dataframe of the association results to plot, or the results output of ```lp.reduce_gaps()```| - |
|**exon_info**| pd.DataFrame: A pandas dataframe of the filtered exon locations, or the exons output of ```lp.reduce_gaps()```| - |
|**title**| string: Title of the plot| - |
|fig_height| int: Figure height (can be adjusted using fig.update_layout()) | 500|
|fig_width | int: Figure height (can be adjusted using fig.update_layout()) | 1400|
|ex_start_col | string: Name of column of exon_info of exome start locations | 'Exon start (bp)'|
|ex_stop_col | string: Name of column of exon_info of exome stop locations| 'Exon stop (bp)'|
|lolli_x | string: Name of column of results df for x axis | 'GENPOS'|
|lolli_y | string: Name of column of results df for y axis | 'LOG10P'|
|lolli_size | string: Name of column of results df for size of bubbles | 'BETA'|
|lolli_col | string:  Name of column of results df for colour of bubbles | 'MASK'|
|lolli_direction | string: Name of column of reults df for direction of bubbles | 'BETA'|
|lollipop_max_size | float: Maximum size of bubbles (Note this need to be a decimal, not an integer) | 5.|
|lollipop_stem_width | float: Width of the lollipop stems | 0.1|

Returns a lolliplot as a plotly figure.

> Note: the exon_info should be filtered to your desired transcript before inputting into ```lp.reduce_gaps()``` or ```lp.lolliplot_raw()```. This can be done in pandas using the ```.loc``` attribute eg: ```df = df.loc[df['Transcript ID'] == TRANSCRIPT_ID, :]```


### lp.example():
Arguments: **None**

This is an example lolliplot script, using available associations on the Exeter server. This function can be modified and run to get lolliplots out quickly and easily. To save the figures, the ```.show()``` parts should be commented out, and the ```write_html()```should be uncommented. The name of the output file should be changed from *lolliplots_example.html*.
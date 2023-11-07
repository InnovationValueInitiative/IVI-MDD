# The IVI-MDD Value Model

As part of the Open-Source Value Project, this repository contains a collection of R scripts that run the [Innovation and Value Initiative's (IVI's)](https://thevalueinitiative.org/) continuous-time individual patient simulation model for major depressive disorder (MDD) (the IVI-MDD model). The model simulates cost and health outcomes associated with a variety of MDD treatments. The model contains four core treatment classes (selective serotonin reuptake inhibitor [SSRIs], serotonin-norepinephrine reuptake inhibitors [SNRIs], atypical antidepressants, and psychotherapy), as well as an additional treatment class (atypical antipsychotics) utilized in combination therapy, and an add-on treatment class (somatic therapies). The options available for each line of therapy are informed by treatment guidelines of the American Psychological Association and American Psychiatric Association, real-world database studies, and input from a Clinical Expert Panel. The model allows for many unique sequences of treatments, each comprised of up to four separate lines of therapy. Each treatment or line of therapy may include only a single treatment class (monotherapy) or select combinations of two distinct treatment classes (e.g., “SSRI + psychotherapy”). In addition, somatic therapy may be included as an “add-on” to fourth-line therapy. The model is intended to allow users to evaluate the benefits and costs associated with various treatment sequencies in adults (age 18-64 years) in the United States newly diagnosed with MDD by a healthcare provider, from multiple perspectives (i.e., private and public payers, employers, people with MDD, and society). Rather than identifying a single set of structural assumptions, the model incorporates the flexibility of including multiple scientifically defensible assumptions, allowing for exploration of structural uncertainty and customization based on user preferences and available data. The intention was that such an approach would help highlight existing data and method gaps to promote conversations across stakeholders and underline areas for future research.

## Documentation

### Technical Report

A copy of the model technical report with complete details on model assumptions can be found [here](LINKREQUIRED).

### General Model Use

The below items should be saved to the user's applicable local working directory to run the model offline:

#### Required packages

The R scripts contain applicable library calls to load required packages. Before running the scripts ensure the appropriate packages are installed on the user's machine. This can be done by running the following in the R console:

```r
# Install packages required to run the IVI-MDD model
install.packages(c('dampack','data.table','doParallel','dplyr','dtplyr','openxlsx','readxl','rjson'))
```

#### ivi_mdd_model.Rproj

The IVI-MDD model scripts are coded such that an appropriate working directory will automatically be set if the scripts are kept in the same directory as the `ivi_mdd_model.Rproj` file. To ensure the scripts are referencing the correct working directory, the end user should open `ivi_mdd_model.Rproj` in an instance of RStudio before opening and running any of the IVI-MDD model scripts. For information on RStudio project files, plus refer to the applicable [RStudio documentation](https://support.posit.co/hc/en-us/articles/200526207-Using-RStudio-Projects).

***Note that RStudio is not required to run the model. The end user may customize the working directories and paths referenced in the appliable IVI-MDD model R scripts for use without the aforementioned RStudio project file.***

#### IVI-MDD model scripts

The R script components of the model are described below:

##### dsa.R

Functions relating to conducting the univariate sensitivity analysis. This script is sourced by the `import_inputs.R` script and simply needs to be placed into the user's working directory; it does not need to be run independently.

##### generate_inputs_json.R

Turns the Microsoft Excel inputs workbook into a JSON file for ingestion of input values into the model. The end user may make changes to inputs within the Excel workbook and run the script to create a custom version of model input values. 

***This script additionally is where the user may set the treatment pathways of interest.*** A user may modify treatment pathways by setting the applicable treatments of interest applied to the `txPaths` list object as well as set the use of fourth line somatic add-on therapy in the `ln4.addonsom.YNs` list object. Note that while five treatment pathways are set by default in this script, the user may increase or decrease the number of treatment pathways as the user sees fit by modifiying the number of pathways entered into the `txPaths` and `ln4.addonsom.YNs` list objects. Treatment abbreviations for use in the `txPaths` list object may be found on lines 40 through 53 of the `import_inputs.R` script; inputs to the `ln4.addonsom.YNs` list object are simply `"Yes"` or `"No"`.

If the user wishes to maintain the default model inputs, this script does not need to be run.

##### import_inputs.R

***Imports the inputs JSON file, runs the model, and runs the univariate sensitivity analysis.*** The end user may modify this script to export model results to either JSON or Excel within the applicable calls of the `f.runModel` and `f.runDSA` functions.

##### ivi_mdd_sim.R

Functions relating to the underlying engine of the model. This script is sourced by the `import_inputs.R` script and simply needs to be placed into the user's working directory; it does not need to be run independently.

##### process_inputs.R

Functions relating to the processing of inputs into the model. This script is sourced by the `import_inputs.R` script and simply needs to be placed into the user's working directory; it does not need to be run independently.

##### summarize_results.R

Functions relating to summarizing and exporting of model results. This script is sourced by the `import_inputs.R` script and simply needs to be placed into the user's working directory; it does not need to be run independently.

#### IVI-MDD model directories

##### inputs

Contains the default Excel inputs workbook and JSON inputs files. The inputs in the Excel workbook may be modified by the user as seen fit and a new JSON inputs file generated using the `generate_inputs_json.R` script (see above section for details).

##### results

The default location for model results exported using the `import_inputs.R` script (see above for details).

## Web Interface

Users may also run the model online using our web interface [here](LINKREQUIRED).
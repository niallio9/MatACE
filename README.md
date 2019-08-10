# MatACE

A collection of tools for reading, writing, exploring, manipulating, and plotting data from the ACE (Atmospheric Chemistry Experiment) satellite. There are two instruments aboard the satellite: ACE-FTS and ACE-MAESTRO.

## Getting Started

The majority (if not all) of the code here is written in MatLAB_R2016b. Matlab is required to run the code.

### Prerequisites

All required tools should be found in within the project folder, or be included with MatLAB. Please let me know if not.

ACE Level 2 data is currently available from https://databace.scisat.ca/level2/. You will need to make an account to access ACE-FTS data. It is recommended to download and use the netCDF version of the data, as opposed to the ASCII version.

It is also recommended to use the ACE-FTS that has flag information with the file. These files are available for v3.6/3.6 of the data, but not currently for v4.0.

```
ACE-FTS data is grouped by gas species. For example, the file for 'ozone' Level 2 data is called: 'ACEFTS_L2_v3p6_O3.nc'.
```

### Installing

The scripts and functions can be used after downloading.

### Examples

To read the above netCDF file for ozone from your current directory:

```
tanstruct = read_ace_ncdata(ACEFTS_L2_v3p6_O3.nc')
```

To read all files from a directory and store them in another directory as matlab structures:

```
read_ace_ncdata_for_mat('all')
```

Note that the above script requires that the input and output directories be edited in the 'Define some things' section at the beginning of the script.


One the data is read into matlab format, it can be used by the other functions.

An example would be to make a global ozone climatology for the year 2010, after including the geolocation (GLC) data:

```
tanstruct = read_ace_ncdata(ACEFTS_L2_v3p6_O3.nc'); % read in the ACE ozone data

glcstruct = read_ace_ncdata_glc(ACEFTS_L2_v3p6_GLC.nc') % read in the geolocation data

tanstruct = merge_ace_glc(tanstruct, glcstruct); clear glcstruct; % merge the GLC data and clear it from the workspace

tanstruct = subset_ace_by_year(tanstruct, 2010); % subset to only data from 2010

climstruct = make_ace_climatology(tanstruct); % make a global climatology of the data
``` 

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

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc






matACE README - NJR 04/2018

These are a group of matlab functions that can be used to read and
manipulate ACE-FTS data. There are also some functions for reading and 
using the ACE DMP information. The ACE data netcdf files and ACE DMP files 
are available online.

The trace species geolocation data for ace occultations are read into
separate matlab structures, and can be merged together afterwards
(basically just adds the latitude and longitude information to the trace
gas structure).

Anyway, below is a list of the current functions (20-04-2018). I reckon I haven't
exhausted error possibilities, so if you notice anything wrong or off, just
let me know. And also if yuou notice that the header info is wrong/unclear/unfinished.


apply_ace_flags.m                          
bin_ace_by_eql.m                           
bin_ace_by_lat.m                           
filter_ace_bad_lat.m                       
filter_ace_bad_lon.m                       
filter_ace_pressure.m                      
get_LST_stats.m                            
get_ace_lst.m                              
get_ace_lst_tangent.m                      
get_ace_occultation_names.m                
include_ace_scaled_apriori.m               
interpolate_ace_dmp_to_altgrid.m           
interpolate_ace_dmp_to_pgrid.m             
interpolate_ace_to_altgrid.m               
interpolate_ace_to_pgrid.m                 
make_ace_NOx_climatology.m                 
make_ace_NOy_climatology.m                 
make_ace_climatology.m                     
make_ace_climatology_3month.m              
make_ace_climatology_3month_by_eql.m       
make_ace_climatology_by_eql.m              
make_ace_climatology_month.m               
make_ace_climatology_month_by_eql.m        
make_ace_climatology_multiple.m            
make_ace_climatology_multiple_by_eql.m     
make_ace_climatology_serialmonth.m         
make_ace_climatology_serialmonth_by_eql.m  
make_ace_empty_eqlbinstruct.m              
make_ace_empty_latbinstruct.m              
make_ace_empty_netcdf_SDI.m                
make_ace_gas_vmrs_with_pratmo.m            
make_ace_nitrogen_ratios_with_pratmo.m     
match_ace_data.m                           
match_ace_data_dmp.m                       
merge_ace_dmp.m                            
merge_ace_glc.m                            
plot_ace_climatology.m                     
plot_ace_climatology_file.m                
plot_ace_climatology_to_pdf.m              
read_ace_dmp_for_mat.m                     
read_ace_dmpv2.m                           
read_ace_ncdata.m                          
read_ace_ncdata_for_mat.m                  
read_ace_ncdata_for_mat_v3p5.m             
read_ace_ncdata_glc.m                      
read_ace_ncdata_v3p5.m                     
reduce_climstruct_data_by_obs_nr.m         
reduce_dmpstruct_by_rowindex.m             
reduce_dmpstruct_data_by_index.m           
reduce_glcstruct_by_rowindex.m             
reduce_tanstruct_by_rowindex.m             
reduce_tanstruct_data_by_index.m           
remove999_ace_dmp.m                        
remove_ace_surface_pressure.m              
scale_ace_with_nitrogen_ratios.m           
scale_ace_with_pratmo.m                    
split_ace_by_lst.m                         
split_ace_by_lst_tangent.m                 
subset_ace_by_3month.m                     
subset_ace_by_date.m                       
subset_ace_by_eql.m                        
subset_ace_by_lat.m                        
subset_ace_by_lst.m                        
subset_ace_by_lst_tangent.m                
subset_ace_by_month.m                      
subset_ace_by_year.m                       
subset_ace_dmp_by_date.m                   
subset_ace_dmp_by_month.m                  
subset_ace_dmp_by_year.m                                                       
write_ace_am_pm.m                          
write_ace_climatology_to_netcdf.m          
write_ace_scaled_apriori.m                 
write_ace_scaled_with_pratmo.m      
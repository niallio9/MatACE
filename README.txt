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
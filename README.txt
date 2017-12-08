matACE v1 README - NJR 11/2017

These are a group of matlab functions that can be used to read and
manipulate ACE-FTS data. There are also some functions for reading and 
using the ACE DMP information. The ACE data netcdf files and ACE DMP files 
are available online.

The trace species geolocation data for ace occultations are read into
separate matlab structures, and can be merged together afterwards
(basically just adds the latitude and longitude information to the trace
gas structure).

Anyway, below is a list of the current functions. I reckon I haven't
exhausted error possibilities, so if you notice anything wrong or off, just
let me know.


apply_ace_flags.m                 read_ace_ncdata_for_mat.m
bin_ace_by_lat.m                  read_ace_ncdata_glc.m
filter_ace_pressure.m             reduce_dmpstruct_by_rowindex.m
get_LST_stats.m                   reduce_dmpstruct_data_by_index.m
get_ace_lst.m                     reduce_glcstruct_by_rowindex.m
get_ace_lst_tangent.m             reduce_tanstruct_by_rowindex.m
get_ace_occultation_names.m       reduce_tanstruct_data_by_index.m
interpolate_ace_dmp_to_altgrid.m  remove999_ace_dmp.m
interpolate_ace_dmp_to_pgrid.m    remove_ace_surface_pressure.m
interpolate_ace_to_altgrid.m      subset_ace_by_date.m
interpolate_ace_to_pgrid.m        subset_ace_by_lat.m
make_ace_climatology.m            subset_ace_by_lst.m
make_ace_empty_latbinstruct.m     subset_ace_by_lst_tangent.m
match_ace_data.m                  subset_ace_by_month.m
match_ace_data_dmp.m              subset_ace_by_year.m
merge_ace_glc.m                   subset_ace_dmp_by_date.m
read_ace_dmp_for_mat.m            subset_ace_dmp_by_month.m
read_ace_dmpv2.m                  subset_ace_dmp_by_year.m
read_ace_ncdata.m
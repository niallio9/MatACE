function [ ] = make_ace_NOy_climatology( )
 %A function to create an ACE climatology for NOx. The function uses .mat
 % ACE climatology files and uses a simple addition of the respective
 % climatologies. The individual files are first filtered by the number of
 % observations used in each climatology. Historically, 5 measurements are
 % required for a grid point in the climtology to be valid.
%
% *INPUT*
%
%
% *OUTPUT*
%           .mat files of the climatology information will be written to
%           'climdir', the same directory as the rest of the climatology
%           files, as defined below by the user.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 03/18

%% define some things
%USER DEFINED
% climdir = '/Users/niall/Dropbox/climatology/nryan/climdata/'; % edit this to your directory that contains the ACE netcdf data
% climdir = '/Volumes/Seagate Backup Plus Drive/ACE/climdata/';
climdir = 'F:\ACE\climdata\';
% climdir = '/net/deluge/pb_1/users/nryan/ACE/climdata/';
if ~isdir(climdir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',climdir)
    error('The directory containing the .mat climatology data couldn''t be found')
end

%STANDARD
gasfolder = {'NO', 'NO2', 'HNO3', 'ClONO2', 'N2O5', 'HNO4' }; % an array with the gas names;
lgas = length(gasfolder);
min_obs = 5; % the minimum number of observations needed to include a climatology point

%% go through the climatology data and pull out whats needed to make the NOx climatology (the NO and NO2 files)

gasdir_no = fullfile(climdir,gasfolder{1}); % make the full path to the folder for NO
gasdir_no2 = fullfile(climdir,gasfolder{2}); % make the full path to the folder for NO
gasdir_hno3 = fullfile(climdir,gasfolder{3}); % make the full path to the folder for NO
gasdir_clono2 = fullfile(climdir,gasfolder{4}); % make the full path to the folder for NO
gasdir_n2o5 = fullfile(climdir,gasfolder{5}); % make the full path to the folder for NO
gasdir_hno4 = fullfile(climdir,gasfolder{6}); % make the full path to the folder for NO
gasdir_noy = fullfile(climdir,'NOy'); % make the full path to the folder for NO
if exist(gasdir_no,'dir') ~=7
    error('there is no folder found called %s. Stopping.', gasdir_no)
end
if exist(gasdir_no2,'dir') ~=7
    error('there is no folder found called %s. Stopping.', gasdir_no2)
end
if exist(gasdir_hno3,'dir') ~=7
    error('there is no folder found called %s. Stopping.', gasdir_hno3)
end
if exist(gasdir_clono2,'dir') ~=7
    error('there is no folder found called %s. Stopping.', gasdir_clono2)
end
if exist(gasdir_n2o5,'dir') ~=7
    error('there is no folder found called %s. Stopping.', gasdir_n2o5)
end
if exist(gasdir_hno4,'dir') ~=7
    error('there is no folder found called %s. Stopping.', gasdir_hno4)
end
if exist(gasdir_noy,'dir') ~=7
    mkdir(gasdir_noy);
end
if exist(fullfile(gasdir_noy,'serial_month'),'dir') ~=7
    mkdir(fullfile(gasdir_noy,'serial_month'));
end
templat_serial = dir(fullfile(gasdir_no,'serial_month','*lat*.mat')); % get the names of the serial-month latitude NO climatology files
climfile_lat_serial_no = {templat_serial.name};
if ~isempty(climfile_lat_serial_no)
    % get the years and months that the NO measurements are for
    for j = 1:length(climfile_lat_serial_no) %
        yearandmonth{j} = climfile_lat_serial_no{j}(end-10:end-4); %#ok<AGROW> % will be small vector
    end
    for j = 1:length(yearandmonth) % go through each file NO file available and see if there is an NO2 file for the same year/month
        %the .mat climatolgy files are of the type
        %"ACEFTS_CLIM_v3_lat_NO2_2013_06.mat"
        no_yearandmonth = sprintf('ACEFTS_CLIM_v3_lat_%s_%s.mat',gasfolder{1},yearandmonth{j});
        climfile_no = fullfile(gasdir_no,'serial_month', no_yearandmonth);
        no2_yearandmonth = sprintf('ACEFTS_CLIM_v3_lat_%s_%s.mat',gasfolder{2},yearandmonth{j});
        climfile_no2 = fullfile(gasdir_no2,'serial_month', no2_yearandmonth); % the name of a climatology file
        hno3_yearandmonth = sprintf('ACEFTS_CLIM_v3_lat_%s_%s.mat',gasfolder{3},yearandmonth{j});
        climfile_hno3 = fullfile(gasdir_hno3,'serial_month', hno3_yearandmonth); % the name of a climatology file
        clono2_yearandmonth = sprintf('ACEFTS_CLIM_v3_lat_%s_%s.mat',gasfolder{4},yearandmonth{j});
        climfile_clono2 = fullfile(gasdir_clono2,'serial_month', clono2_yearandmonth); % the name of a climatology file
        n2o5_yearandmonth = sprintf('ACEFTS_CLIM_v3_lat_%s_%s.mat',gasfolder{5},yearandmonth{j});
        climfile_n2o5 = fullfile(gasdir_n2o5,'serial_month', n2o5_yearandmonth); % the name of a climatology file
        hno4_yearandmonth = sprintf('ACEFTS_CLIM_v3_lat_%s_%s.mat',gasfolder{6},yearandmonth{j});
        climfile_hno4 = fullfile(gasdir_hno4,'serial_month', hno4_yearandmonth); % the name of a climatology file
        if exist(climfile_no2,'file') == 2 && exist(climfile_hno3,'file') == 2 && exist(climfile_clono2,'file') == 2 && exist(climfile_n2o5,'file') == 2 && exist(climfile_hno4,'file') ==2 % check if each file exists for that month
            load(climfile_no2); % loads a variable called climstruct
            climstruct_no2 = climstruct;
            load(climfile_no); % loads a variable called climstruct
            climstruct_no = climstruct;
            load(climfile_hno3); % loads a variable called climstruct
            climstruct_hno3 = climstruct;
            load(climfile_clono2); % loads a variable called climstruct
            climstruct_clono2 = climstruct;
            load(climfile_n2o5); % loads a variable called climstruct
            climstruct_n2o5 = climstruct;
            load(climfile_hno4); % loads a variable called climstruct
            climstruct_hno4 = climstruct;
            clear climstruct
            % reduce the climstruct to only include climatology points that are created by more then 5 observations
            climstruct_no = reduce_climstruct_data_by_obs_nr(climstruct_no, min_obs);
            climstruct_no2 = reduce_climstruct_data_by_obs_nr(climstruct_no2, min_obs);
            climstruct_hno3 = reduce_climstruct_data_by_obs_nr(climstruct_hno3, min_obs);
            climstruct_clono2 = reduce_climstruct_data_by_obs_nr(climstruct_clono2, min_obs);
            climstruct_n2o5 = reduce_climstruct_data_by_obs_nr(climstruct_n2o5, min_obs);
            climstruct_hno4 = reduce_climstruct_data_by_obs_nr(climstruct_hno4, min_obs);
            %get the data you need and make the NOx climatology 
            climstruct_noy.source_file = strcat(climstruct_no.source_file,', ', climstruct_no2.source_file,', ', climstruct_hno3.source_file,', ', climstruct_clono2.source_file,', ', climstruct_n2o5.source_file,', ', climstruct_hno4.source_file);
            climstruct_noy.beta_angle_mean = (climstruct_no.beta_angle_mean + climstruct_no2.beta_angle_mean + climstruct_hno3.beta_angle_mean + climstruct_clono2.beta_angle_mean + climstruct_n2o5.beta_angle_mean)/lgas + climstruct_hno4.beta_angle_mean;
            climstruct_noy.doy_mean = (climstruct_no.doy_mean + climstruct_no2.doy_mean + climstruct_hno3.doy_mean + climstruct_clono2.doy_mean + climstruct_n2o5.doy_mean + climstruct_hno4.doy_mean)/lgas;
            climstruct_noy.start_date = datestr(datenum(climstruct_no.start_date) + datenum(climstruct_no2.start_date) + datenum(climstruct_hno3.start_date) + datenum(climstruct_clono2.start_date) + datenum(climstruct_n2o5.start_date) + datenum(climstruct_hno4.start_date))/lgas;
            climstruct_noy.end_date = datestr(datenum(climstruct_no.end_date) + datenum(climstruct_no2.end_date) + datenum(climstruct_hno3.end_date) + datenum(climstruct_clono2.end_date) + datenum(climstruct_n2o5.end_date) + datenum(climstruct_hno4.end_date))/lgas;
            climstruct_noy.gas = 'NOy';
            climstruct_noy.lon_tangent_mean = (climstruct_no.lon_tangent_mean + climstruct_no2.lon_tangent_mean + climstruct_hno3.lon_tangent_mean + climstruct_clono2.lon_tangent_mean + climstruct_n2o5.lon_tangent_mean + climstruct_hno4.lon_tangent_mean)/lgas;
            climstruct_noy.lat_tangent_mean = (climstruct_no.lat_tangent_mean + climstruct_no2.lat_tangent_mean + climstruct_hno3.lat_tangent_mean + climstruct_clono2.lat_tangent_mean + climstruct_n2o5.lat_tangent_mean + climstruct_hno4.lat_tangent_mean)/lgas;
            climstruct_noy.altitude_km_mean = (climstruct_no.altitude_km_mean + climstruct_no2.altitude_km_mean + climstruct_hno3.altitude_km_mean + climstruct_clono2.altitude_km_mean + climstruct_n2o5.altitude_km_mean + climstruct_hno4.altitude_km_mean)/lgas;
            climstruct_noy.pressure_hPa = climstruct_no.pressure_hPa; %this should be the same for all gases
            climstruct_noy.vmr_zonal = (climstruct_no.vmr_zonal + climstruct_no2.vmr_zonal + climstruct_hno3.vmr_zonal + climstruct_clono2.vmr_zonal + climstruct_n2o5.vmr_zonal + climstruct_hno4.vmr_zonal)/lgas;
            climstruct_noy.vmr_zonal_var = climstruct_no.vmr_zonal_var + climstruct_no2.vmr_zonal_var + climstruct_hno3.vmr_zonal_var + climstruct_clono2.vmr_zonal_var + climstruct_n2o5.vmr_zonal_var + climstruct_hno4.vmr_zonal_var; % errors added in quadrature is the same as adding variances
            climstruct_noy.vmr_zonal_error = sqrt(climstruct_no.vmr_zonal_error.^2 + climstruct_no2.vmr_zonal_error.^2 + climstruct_hno3.vmr_zonal_error.^2 + climstruct_clono2.vmr_zonal_error.^2 + climstruct_n2o5.vmr_zonal_error.^2 + climstruct_hno4.vmr_zonal_error.^2); % add in quadrature
            climstruct_noy.vmr_zonal_standard_error = sqrt(climstruct_no.vmr_zonal_standard_error.^2 + climstruct_no2.vmr_zonal_error.^2 + climstruct_hno3.vmr_zonal_error.^2 + climstruct_clono2.vmr_zonal_error.^2 + climstruct_n2o5.vmr_zonal_error.^2 + climstruct_hno4.vmr_zonal_error.^2);
            climstruct_noy.lat = climstruct_no.lat; %this should be the same for all gases
            climstruct_noy.lat_bounds = climstruct_no.lat_bounds; %this should be the same for all gases
            climstruct_noy.lon_mean = (climstruct_no.lon_mean + climstruct_no2.lon_mean + climstruct_hno3.lon_mean + climstruct_clono2.lon_mean + climstruct_n2o5.lon_mean + climstruct_hno4.lon_mean)/lgas;
            climstruct_noy.obs_count = (climstruct_no.obs_count + climstruct_no2.obs_count + climstruct_hno3.obs_count + climstruct_clono2.obs_count + climstruct_n2o5.obs_count + climstruct_hno4.obs_count)/lgas;
            climstruct_noy.obs_location = cell(size(climstruct_no.obs_location)); % can't really average these
            climstruct_noy.lst_median = (climstruct_no.lst_median + climstruct_no2.lst_median + climstruct_hno3.lst_median + climstruct_clono2.lst_median + climstruct_n2o5.lst_median + climstruct_hno4.lst_median)/lgas;
            climstruct_noy.lst_mean = (climstruct_no.lst_mean + climstruct_no2.lst_mean + climstruct_hno3.lst_mean + climstruct_clono2.lst_mean + climstruct_n2o5.lst_mean + climstruct_hno4.lst_mean)/lgas;
            climstruct_noy.lst_max = (climstruct_no.lst_max + climstruct_no2.lst_max + climstruct_hno3.lst_max + climstruct_clono2.lst_max + climstruct_n2o5.lst_max + climstruct_hno4.lst_max)/lgas;
            climstruct_noy.lst_min = (climstruct_no.lst_min + climstruct_no2.lst_min + climstruct_hno3.lst_min + climstruct_clono2.lst_min + climstruct_n2o5.lst_min + climstruct_hno4.lst_min)/lgas;
            climstruct_noy.lst_var = climstruct_no.lst_var + climstruct_no2.lst_var + climstruct_hno3.lst_var + climstruct_clono2.lst_var + climstruct_n2o5.lst_var + climstruct_hno4.lst_var;
            climstruct_noy.climatology_type = climstruct_no.climatology_type;
            climstruct_noy.time = climstruct_no.time;
            % save the file
            climstruct = climstruct_noy;
            savedest = fullfile(gasdir_noy,'serial_month',sprintf('ACEFTS_CLIM_v3_lat_NOy_%s.mat', yearandmonth{j}));
            fprintf('\nSaving NOy %s climatology to %s\n', yearandmonth{j}, savedest);
            save(savedest,'climstruct');
        end
    end
end
disp('All done :)')
%
end


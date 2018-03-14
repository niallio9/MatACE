function [ ] = make_ace_NOx_climatology( )
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
climdir = '/Volumes/Seagate Backup Plus Drive/ACE/climdata/';
% climdir = 'F:\ACE\climdata\';
% climdir = '/net/deluge/pb_1/users/nryan/ACE/climdata/';
if ~isdir(climdir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',climdir)
    error('The directory containing the .mat climatology data couldn''t be found')
end

%STANDARD
gasfolder = {'NO', 'NO2'}; % an array with the gas names;
min_obs = 5; % the minimum number of observations needed to include a climatology point

%% go through the climatology data and pull out whats needed to make the NOx climatology (the NO and NO2 files)

gasdir_no = fullfile(climdir,gasfolder{1}); % make the full path to the folder for NO
gasdir_no2 = fullfile(climdir,gasfolder{2}); % make the full path to the folder for NO
gasdir_nox = fullfile(climdir,'NOx'); % make the full path to the folder for NO
if exist(gasdir_no,'dir') ~=7
    error('there is no folder found called %s. Stopping.', gasdir_no)
end
if exist(gasdir_no2,'dir') ~=7
    error('there is no folder found called %s. Stopping.', gasdir_no2)
end
if exist(gasdir_nox,'dir') ~=7
    mkdir(gasdir_nox);
end
if exist(fullfile(gasdir_nox,'serial_month'),'dir') ~=7
    mkdir(fullfile(gasdir_nox,'serial_month'));
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
        if exist(climfile_no2,'file') == 2 % check if the file exists for that month
            load(climfile_no2); % loads a variable called climstruct
            climstruct_no2 = climstruct;
            load(climfile_no); % loads a variable called climstruct
            climstruct_no = climstruct;
            clear climstruct
            % reduce the climstruct to only include climatology points that are created by more then 5 observations
            climstruct_no = reduce_climstruct_data_by_obs_nr(climstruct_no, min_obs);
            climstruct_no2 = reduce_climstruct_data_by_obs_nr(climstruct_no2, min_obs);
            %get the data you need and make the NOx climatology 
            climstruct_nox.source_file = strcat(climstruct_no.source_file,', ', climstruct_no2.source_file);
            climstruct_nox.beta_angle_mean = (climstruct_no.beta_angle_mean + climstruct_no2.beta_angle_mean)/2;
            climstruct_nox.doy_mean = (climstruct_no.doy_mean + climstruct_no2.doy_mean)/2;
            climstruct_nox.start_date = datestr(datenum(climstruct_no.start_date) + datenum(climstruct_no2.start_date))/2;
            climstruct_nox.end_date = datestr(datenum(climstruct_no.end_date) + datenum(climstruct_no2.end_date))/2;
            climstruct_nox.gas = 'NOx';
            climstruct_nox.lon_tangent_mean = (climstruct_no.lon_tangent_mean + climstruct_no2.lon_tangent_mean)/2;
            climstruct_nox.lat_tangent_mean = (climstruct_no.lat_tangent_mean + climstruct_no2.lat_tangent_mean)/2;
            climstruct_nox.altitude_km_mean = (climstruct_no.altitude_km_mean + climstruct_no2.altitude_km_mean)/2;
            climstruct_nox.pressure_hPa = climstruct_no.pressure_hPa; %this should be the same for all gases
            climstruct_nox.vmr_zonal = (climstruct_no.vmr_zonal + climstruct_no2.vmr_zonal)/2;
            climstruct_nox.vmr_zonal_var = climstruct_no.vmr_zonal_var + climstruct_no2.vmr_zonal_var; % errors added in quadrature is the same as adding variances
            climstruct_nox.vmr_zonal_error = sqrt(climstruct_no.vmr_zonal_error.^2 + climstruct_no2.vmr_zonal_error.^2); % add in quadrature
            climstruct_nox.vmr_zonal_standard_error = sqrt(climstruct_no.vmr_zonal_standard_error.^2 + climstruct_no2.vmr_zonal_error.^2);
            climstruct_nox.lat = climstruct_no.lat; %this should be the same for all gases
            climstruct_nox.lat_bounds = climstruct_no.lat_bounds; %this should be the same for all gases
            climstruct_nox.lon_mean = (climstruct_no.lon_mean + climstruct_no2.lon_mean)/2;
            climstruct_nox.obs_count = (climstruct_no.obs_count + climstruct_no2.obs_count)/2;
            climstruct_nox.obs_location = cell(size(climstruct_no.obs_location)); % can't really average these
            climstruct_nox.lst_median = (climstruct_no.lst_median + climstruct_no2.lst_median)/2;
            climstruct_nox.lst_mean = (climstruct_no.lst_mean + climstruct_no2.lst_mean)/2;
            climstruct_nox.lst_max = (climstruct_no.lst_max + climstruct_no2.lst_max)/2;
            climstruct_nox.lst_min = (climstruct_no.lst_min + climstruct_no2.lst_min)/2;
            climstruct_nox.lst_var = climstruct_no.lst_var + climstruct_no2.lst_var;
            climstruct_nox.climatology_type = climstruct_no.climatology_type;
            climstruct_nox.time = climstruct_no.time;
            % save the file
            climstruct = climstruct_nox;
            savedest = fullfile(gasdir_nox,'serial_month',sprintf('ACEFTS_CLIM_v3_lat_NOx_%s.mat', yearandmonth{j}));
            fprintf('\nSaving NOx %s climatology to %s\n', yearandmonth{j}, savedest);
            save(savedest,'climstruct');
        end
    end
end
disp('All done :)')
%
end


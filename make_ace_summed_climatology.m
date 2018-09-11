function [ ] = make_ace_summed_climatology( gasout_name, varargin )
 %A function to create an ACE climatology for a family of gases. The
 % function uses .mat ACE climatology files and uses a simple addition of
 % the respective climatologies.
%
% *INPUT*
%           gasout_name: STRING - the name of the family of gases. This is
%           used to name the output file, etc..
%
%           varargin: STRING - the name(s) of the gases that make up the
%           family of gases. You can input the same gas more than once if
%           it makes a multiple contribution to a family.
%           An input for NOx (NO + NO2) would be: 'NO', 'NO2'.
%
%
% *OUTPUT*
%           .mat files of the climatology information will be written to
%           'climdir', the same directory as the rest of the climatology
%           files, as defined below by the user.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 05/18

%% define some things
%USER DEFINED
% climdir = '/Users/niall/Dropbox/climatology/nryan/climdata/'; % edit this to your directory that contains the ACE netcdf data
% climdir = '/Volumes/Seagate Backup Plus Drive/ACE/climdata/';
% climdir = 'F:\ACE\climdata\';
climdir = 'C:\Users\ryann\ACE\climdata_testing\';
% climdir = '/net/deluge/pb_1/users/nryan/ACE/climdata/';
if ~isdir(climdir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',climdir)
    error('The directory containing the .mat climatology data couldn''t be found')
end

clim_type = ''; % this is where you choose what type of climatology to use: '_am', '_pm', '_sap', '_sap_am', '_sap_pm', or ''.

%STANDARD
gasnames = varargin; % array of cells with gas names
lgases = length(gasnames); % number of gases to read
gasfolders = cell(lgases,1);
gasdirs = cell(lgases,1);
gas_yearandmonth = cell(lgases,1);
climfiles = cell(lgases,1);

for n = 1:lgases
    %get a cell of folder names and paths to the folders
    gasfolders{n} = strcat(gasnames{n}, clim_type);
    gasdirs{n} = fullfile(climdir,gasfolders{n});
    if exist(gasdirs{n},'dir') ~= 7
        error('there is no folder found called %s. Stopping.', gasdirs{n})
    end
    % make an output folder if it doesn't exist
    gasdir_out = fullfile(climdir,gasout_name);
    if exist(gasdir_out,'dir') ~=7
        mkdir(gasdir_out);
    end
    if exist(fullfile(gasdir_out,'serial_month'),'dir') ~=7
        mkdir(fullfile(gasdir_out,'serial_month'));
    end
end
min_obs = 5; %#ok<NASGU> % the minimum number of observations needed to include a climatology point
templat_serial = dir(fullfile(gasdirs{n},'serial_month','*lat*.mat')); %  get the names of the serial-month latitude climatology files, for the forst gas in the list
climfile_lat_serial = {templat_serial.name};
if ~isempty(climfile_lat_serial)
    % get the years and months that the measurements are for
    for j = 1:length(climfile_lat_serial) %
        yearandmonth{j} = climfile_lat_serial{j}(end-10:end-4); %#ok<AGROW> % will be small vector
    end
    for j = 1:length(yearandmonth) %
        %the .mat climatolgy files are of the type
        %"ACEFTS_CLIM_v3_lat_NO2_2013_06.mat"
        do_yearandmonth = 1;
        % go through each file available and see if there is an equivalent file for the same year/month for all gases
        for n = 1:lgases
            gas_yearandmonth{n} = sprintf('ACEFTS_CLIM_v3_lat_%s_%s.mat', gasfolders{n}, yearandmonth{j});
            climfiles{n} = fullfile(gasdirs{n},'serial_month', gas_yearandmonth{n}); % the name of a climatology file
            if ~exist(climfiles{n},'file')
                do_yearandmonth = 0;
            end
        end
        if do_yearandmonth == 1 % if all the required files are available
            % create n variables with the data
            for n = 1:lgases
                load(climfiles{n}); % loads a variable called climstruct
                loadgas = sprintf('climstruct_%i = climstruct;', n); % use a numbered system here so you can use a gas more than once
                eval(loadgas);
            end
            clear climstruct
            % reduce the climstruct to only include climatology points that are created by more then 5 observations
            for n = 1:lgases
                reducegas = sprintf('climstruct_%i = reduce_climstruct_data_by_obs_nr(climstruct_%i, min_obs);', n, n);
                eval(reducegas);
            end
            %reset the fields of the output structure
            gasout.source_file = [];
            gasout.beta_angle_mean = 0;
            gasout.doy_mean = 0;
            gasout.start_date = 0;
            gasout.end_date = 0;
%             gasout.gas = [];
            gasout.lon_tangent_mean = 0;
            gasout.lat_tangent_mean = 0;
            gasout.altitude_km_mean = 0;
%             gasout.pressure_hPa = []; %this should be the same for all gases
            gasout.vmr_zonal = 0;
            gasout.vmr_zonal_var = 0;
            gasout.vmr_zonal_error = 0;
            gasout.vmr_zonal_standard_error = 0;
            gasout.lat = 0;
            gasout.lat_bounds = 0;
            gasout.lon_mean = 0;
            gasout.obs_count = 0;
%             gasout.obs_location = {};
            gasout.lst_median = 0;
            gasout.lst_mean = 0;
            gasout.lst_max = 0;
            gasout.lst_min = 0;
            gasout.lst_var = 0;
%             gasout.climatology_type = [];
%             gasout.time = [];
            %get the data you need and make the NOx climatology
            for n = 1:lgases
%                 n
                eval(sprintf('source_file = climstruct_%i.source_file;', n)); % info for gas n.
                gasout.source_file = strcat(gasout.source_file,',', source_file); % append info for gas n
                eval(sprintf('beta_angle_mean = climstruct_%i.beta_angle_mean;', n)); % info for gas n.
                gasout.beta_angle_mean = gasout.beta_angle_mean + beta_angle_mean; % add info for gas n. will devide by n later
                eval(sprintf('doy_mean = climstruct_%i.doy_mean;', n)); % info for gas n.
                gasout.doy_mean = gasout.doy_mean + doy_mean;
                eval(sprintf('start_date = climstruct_%i.start_date;', n)); % info for gas n.
                gasout.start_date = gasout.start_date + datenum(start_date); % will datestr this later
                eval(sprintf('end_date = climstruct_%i.end_date;', n)); % info for gas n.
                gasout.end_date = gasout.end_date + datenum(end_date); % will datestr this later
                gasout.gas = gasout_name;
                eval(sprintf('lon_tangent_mean = climstruct_%i.lon_tangent_mean;', n)); % info for gas n.
                gasout.lon_tangent_mean = gasout.lon_tangent_mean + lon_tangent_mean; % will devide by n later
                eval(sprintf('lat_tangent_mean = climstruct_%i.lat_tangent_mean;', n)); % info for gas n.
                gasout.lat_tangent_mean = gasout.lat_tangent_mean + lat_tangent_mean; % will devide by n later
                eval(sprintf('altitude_km_mean = climstruct_%i.altitude_km_mean;', n)); % info for gas n.
                gasout.altitude_km_mean = gasout.altitude_km_mean + altitude_km_mean; % will devide by n later
                gasout.pressure_hPa = climstruct_1.pressure_hPa; % this should be the same for all gases
                eval(sprintf('vmr_zonal = climstruct_%i.vmr_zonal;', n)); % info for gas n.
                gasout.vmr_zonal = gasout.vmr_zonal + vmr_zonal;
                eval(sprintf('vmr_zonal_var = climstruct_%i.vmr_zonal_var;', n)); % info for gas n.
                gasout.vmr_zonal_var = gasout.vmr_zonal_var + vmr_zonal_var;
                eval(sprintf('vmr_zonal_error = climstruct_%i.vmr_zonal_error;', n)); % info for gas n.
                gasout.vmr_zonal_error = gasout.vmr_zonal_error + vmr_zonal_error.^2; % will square root later
                eval(sprintf('vmr_zonal_standard_error = climstruct_%i.vmr_zonal_standard_error;', n)); % info for gas n.
                gasout.vmr_zonal_standard_error = gasout.vmr_zonal_standard_error + vmr_zonal_standard_error.^2; % will square root later
                gasout.lat = climstruct_1.lat; % this should be the same for all gases
                gasout.lat_bounds = climstruct_1.lat_bounds; % this should be the same for all gases
                eval(sprintf('lon_mean = climstruct_%i.lon_mean;', n)); % info for gas n.
                gasout.lon_mean = gasout.lon_mean + lon_mean; % will devide by n later
                eval(sprintf('obs_count = climstruct_%i.obs_count;', n)); % info for gas n.
                gasout.obs_count = gasout.obs_count + obs_count; % will devide by n later
%                 eval(sprintf('obs_location = climstruct_%i.obs_location;', n)); % info for gas n.
%                 gasout.obs_location = gasout.obs_location + obs_location; % will devide by n later
                gasout.obs_location = cell(size(climstruct_1.obs_location)); % can't really average these
                eval(sprintf('lst_median = climstruct_%i.lst_median;', n)); % info for gas n.
                gasout.lst_median = gasout.lst_median + lst_median; % will devide by n later
                eval(sprintf('lst_mean = climstruct_%i.lst_mean;', n)); % info for gas n.
                gasout.lst_mean = gasout.lst_mean + lst_mean; % will devide by n later
                eval(sprintf('lst_max = climstruct_%i.lst_max;', n)); % info for gas n.
                gasout.lst_max = gasout.lst_max + lst_max; % will devide by n later
                eval(sprintf('lst_min = climstruct_%i.lst_min;', n)); % info for gas n.
                gasout.lst_min = gasout.lst_min + lst_min; % will devide by n later
                eval(sprintf('lst_var = climstruct_%i.lst_var;', n)); % info for gas n.
                gasout.lst_var = gasout.lst_var + lst_var;
                gasout.climatology_type = climstruct_1.climatology_type; % this should be the same for all gases
                gasout.time = climstruct_1.time; % this should be the same for all gases
            end
            % make the means and sums in quadrature that are required for some
            % of the fields.
            gasout.beta_angle_mean = gasout.beta_angle_mean/lgases;
            gasout.doy_mean = gasout.doy_mean/lgases;
            gasout.start_date = gasout.start_date/lgases;
            gasout.end_date = gasout.end_date/lgases;
            gasout.lon_tangent_mean = gasout.lon_tangent_mean/lgases;
            gasout.lat_tangent_mean = gasout.lat_tangent_mean/lgases;
            gasout.altitude_km_mean = gasout.altitude_km_mean/lgases;
            gasout.vmr_zonal_error = sqrt(gasout.vmr_zonal_error);
            gasout.vmr_zonal_standard_error = sqrt(gasout.vmr_zonal_standard_error);
            gasout.lon_mean = gasout.lon_mean/lgases;
            gasout.obs_count = gasout.obs_count/lgases;
%             gasout.obs_location = gasout.obs_location/lgases;
            gasout.lst_median = gasout.lst_median/lgases;
            gasout.lst_mean = gasout.lst_mean/lgases;
            gasout.lst_max = gasout.lst_max/lgases;
            gasout.lst_min = gasout.lst_min/lgases;
            % save the file
            climstruct = gasout; %#ok<NASGU>
            savedest = fullfile(gasdir_out,'serial_month',sprintf('ACEFTS_CLIM_v3_lat_%s_%s.mat', gasout_name, yearandmonth{j}));
            fprintf('\nSaving %s %s climatology to %s\n', gasdir_out, yearandmonth{j}, savedest);
            save(savedest,'climstruct');
        else
            fprintf('\nThe files for each gas aren''t available for %s', yearandmonth{j})
        end
        
    end
end
   
disp('All done :)')
%
end


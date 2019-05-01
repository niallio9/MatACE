function [ ] = write_ace_climatology_to_netcdf_SDI( varargin )
 %A function to create an ACE climatology netCDF file from the information
 %contained in the ACE climstruct .mat file.
%
% *INPUT*
%
%           varargin: STRING - the name of the gas for which you want to write a
%           netcdf file from the .mat climatology file.
%           The input may be a single gas (e.g., 'O3'), or multiple
%           gases (e.g., 'O3', 'ClO').
%           To read all .mat files in the climatology data directory, the
%           input is 'all'.
%
% *OUTPUT*
%           .nc files of the climatology information will be written to
%           'ncdir', as defined below by the user.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 03/18

%% define some things
%USER DEFINED
% climdir = '/Users/niall/Dropbox/climatology/nryan/climdata/'; % edit this to your directory that contains the ACE netcdf data
% climdir = '/Volumes/Seagate Backup Plus Drive/ACE/climdata/';
climdir = 'C:\Users\ryann\ACE\climdata\scaled_with_all_data';
% climdir = 'F:\ACE\climdata\';
% climdir = '/net/deluge/pb_1/users/nryan/ACE/climdata/';
if ~isdir(climdir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',climdir)
    error('The directory containing the .mat climatology data couldn''t be found')
end
% ncdir = '/Users/niall/Dropbox/climatology/nryan/climdata_netcdf/'; % edit this to your output directory
% ncdir = '/Volumes/Seagate Backup Plus Drive/ACE/climdata_netcdf/';
ncdir = 'C:\Users\ryann\ACE\climdata_netcdf';
% ncdir = 'F:\ACE\climdata_netcdf\';
% ncdir = '/net/deluge/pb_1/users/nryan/ACE/climdata_netcdf/';

%STANDARD
gas_in = string(varargin); % an array with the gas names;
lin = length(gas_in);
min_obs = 5; % the minimum number of observations needed to include a climatology point
% llev = 48;
llev = 28; levind = 6:33;
llat = 36;
ltime = 12;
% vmr_out = -999*ones(llat,llev, ltime);
% vmr_std_out = -999*ones(llat,llev, ltime);
% vmr_obs_out = -999*ones(llat,llev, ltime);
% lst_mean = -999*ones(llat, ltime);
% lst_max = -999*ones(llat, ltime);
% lst_min = -999*ones(llat, ltime);
% ave_dom = -999*ones(llat, ltime);
% ave_lat = -999*ones(llat, ltime);
% lat = -999*ones(llat,1);
% plev = -999*ones(llev,1);
% time = -999*ones(ltime,1);
sdate1950 = datenum(1950,1,1); % the serial date of 1950-01-01 00:00:00


%% go through the climatology data and pull out whats needed to make the netcdf file
% for the serail month data you'll need to read all the data for one year
% and put it into a single netcdf file
if isdir(ncdir)
    temp = dir(climdir); % get the names of the folders (named by gas) that contain the data
    temp = {temp.name};
    % remove folders that start with '.'
    gasfolder = cell(150,1); % will remove the empty cells later
    for i = 1:length(temp)
        if ~strcmp(temp{i}(1), '.')
            gasfolder{i} = temp{i};
        end
    end
    gasfolder = gasfolder(~cellfun('isempty',gasfolder)); % remove the empty cells
    gasfolder = string(gasfolder); % changed this now to make choosing a gas with main function input easier.
    % reduce the gasfolder names to only include those mentioned in the
    % main function input. Nothing is done if the input is 'all'.
    if lin == 1
        if ~strcmp(gas_in,'all')
            gasfolder = intersect(gasfolder, gas_in);
        end
    else
        gasfolder = intersect(gasfolder, gas_in); 
    end
    
    for i = 1:length(gasfolder) % loop through the gases 
        gasdir_i = fullfile(climdir,gasfolder{i}); % make the full path to the folder for a gas
        if exist(gasdir_i,'dir') ~=7
            error('there is no folder found called %s. Stopping.', gasdir_i)
        end
        templat_serial = dir(fullfile(gasdir_i,'serial_month','*lat*.mat')); % get the names of the serial-month latitude climatology files
        climfile_lat_serial = {templat_serial.name};
        if ~isempty(climfile_lat_serial)
            % get the years that the measurements are for
            for j = 1:length(climfile_lat_serial) %
                data_years(j) = str2double(climfile_lat_serial{j}(end-10:end-7)); %#ok<AGROW> % will be small vector
            end
            data_years_unique = unique(data_years);
            for j = 1:length(data_years_unique) % make a netcdf file for each year of data
                % set/reset the variables as they are filled by
                % index for each year
                vmr_out = -999*ones(llat,llev, ltime);
                vmr_std_out = -999*ones(llat,llev, ltime);
                vmr_obs_out = -999*ones(llat,llev, ltime);
                lst_mean = -999*ones(llat, ltime);
                lst_max = -999*ones(llat, ltime);
                lst_min = -999*ones(llat, ltime);
                ave_dom = -999*ones(llat, ltime);
                ave_lat = -999*ones(llat, ltime);
                lat = -999*ones(llat,1);
                plev = -999*ones(llev,1);
                time = -999*ones(ltime,1);
                % get the data for each month of the given year
                for k = 1:12
                    %the .mat climatolgy files are of the type
                    %"ACEFTS_CLIM_v3_eql_NO2_2013_06.mat"
                    climfile_yearandmonth = sprintf('ACEFTS_CLIM_v3_lat_%s_%i_%02.0f.mat',gasfolder{i},data_years_unique(j),k);
                    climfile_k = fullfile(gasdir_i,'serial_month', climfile_yearandmonth); % the name of a climatology file
                    if exist(climfile_k,'file') == 2 % check if the file exists for that month
                        load(climfile_k); % loads a variable called climstruct
                        % reduce the climstruct to only include climatology points that are created by more then 5 observations
                        if length(gasfolder{i}) > 2
                            if ~strcmp(gasfolder{i}(1:3), 'NOy') || ~strcmp(gasfolder{i}(1:3), 'NOx')
                            climstruct = reduce_climstruct_data_by_obs_nr(climstruct, min_obs);
                            end
                        end
                        %get the data you need for the netcdf file
                        vmr_out(:,:,k) = climstruct.vmr_zonal(levind,:)'; %36x48 array
                        vmr_std_out(:,:,k) = sqrt(climstruct.vmr_zonal_var(levind,:)'); % want the standard deviation here
                        vmr_obs_out(:,:,k) = climstruct.obs_count(levind,:)';
                        lst_mean(:,k) = climstruct.lst_mean;
                        lst_max(:,k) = climstruct.lst_max;
                        lst_min(:,k) = climstruct.lst_min;
                        repyear = repmat(climstruct.time(1),1,length(climstruct.doy_mean));
                        ave_dom(:,k) = day(doy2date(climstruct.doy_mean, repyear));
                        ave_lat(:,k) = climstruct.lat_tangent_mean;

%                         lat = climstruct.lat;
                        
                        time(k) = datenum(data_years_unique(j), k, 15) - sdate1950;
%                         vmr_out(k,:,:) = climstruct.vmr_zonal; %48x36 array
%                         vmr_out(isnan(vmr_out)) = -999; % change the nans to -999 values. ***should make this more efficient here***
%                         vmr_std_out(k,:,:) = sqrt(climstruct.vmr_zonal_var); % want the standard deviation here
%                         vmr_std_out(isnan(vmr_std_out)) = -999;
%                         vmr_obs_out(k,:,:) = climstruct.obs_count;
%                         vmr_obs_out(isnan(vmr_obs_out)) = -999;
%                         lst_mean(k,:) = climstruct.lst_mean;
%                         lst_mean(isnan(lst_mean)) = -999;
%                         lst_max(k,:) = climstruct.lst_max;
%                         lst_max(isnan(lst_max)) = -999;
%                         lst_min(k,:) = climstruct.lst_min;
%                         lst_min(isnan(lst_min)) = -999;
%                         repyear = repmat(climstruct.time(1),1,length(climstruct.doy_mean));
%                         ave_dom(k,:) = day(doy2date(climstruct.doy_mean, repyear));
%                         ave_dom(isnan(ave_dom)) = -999;
%                         ave_lat(k,:) = climstruct.lat_tangent_mean;
%                         ave_lat(isnan(ave_lat)) = -999;
%                         lat(:) = climstruct.lat;
%                         plev(:) = climstruct.pressure_hPa;
%                         time(k) = datenum(data_years_unique(j), k, 15) - sdate1950;
                    end
                end
                vmr_out = flip(vmr_out,1); % flip the latitude dimension to be from north to south
                vmr_out(isnan(vmr_out)) = -999; % change the nans to -999 values. ***should make this more efficient here***
                vmr_std_out = flip(vmr_std_out,1); % flip the latitude dimension to be from north to south
                vmr_std_out(isnan(vmr_std_out)) = -999;
                vmr_obs_out = flip(vmr_obs_out,1); % flip the latitude dimension to be from north to south
                vmr_obs_out(isnan(vmr_obs_out)) = -999;
                lst_mean = flip(lst_mean,1); % flip the latitude dimension to be from north to south
                lst_mean(isnan(lst_mean)) = -999;
                lst_max = flip(lst_max,1); % flip the latitude dimension to be from north to south
                lst_max(isnan(lst_max)) = -999;
                lst_min = flip(lst_min,1); % flip the latitude dimension to be from north to south
                lst_min(isnan(lst_min)) = -999;
                ave_dom = flip(ave_dom,1); % flip the latitude dimension to be from north to south
                ave_dom(isnan(ave_dom)) = -999;
                ave_lat = flip(ave_lat,1); % flip the latitude dimension to be from north to south
                ave_lat(isnan(ave_lat)) = -999;
                
                lat(:) = flip(climstruct.lat); % flip the latitude dimension to be from north to south
                plev(:) = climstruct.pressure_hPa(levind,:);
                
                
                % get the name of the gas for the output file
                if length(gasfolder{i}) > 10 && strcmp(gasfolder{i}(end-8:end-6), 'sap') % for '_sap_s10am' and '_sap_s10pm', etc.
                    gasfolder_short = gasfolder{i}(1:end-10);
                elseif length(gasfolder{i}) > 7 && strcmp(gasfolder{i}(end-5:end), 'sap_am') % for '_sap_am' and '_sap_pm'
                    gasfolder_short = gasfolder{i}(1:end-7);
                elseif length(gasfolder{i}) > 7 && strcmp(gasfolder{i}(end-5:end), 'sap_pm') % for '_sap_am' and '_sap_pm'
                    gasfolder_short = gasfolder{i}(1:end-7);
                elseif length(gasfolder{i}) > 4 && strcmp(gasfolder{i}(end-2:end), 'sap') % for '_sap'
                    gasfolder_short = gasfolder{i}(1:end-4);
                elseif length(gasfolder{i}) > 3 && strcmp(gasfolder{i}(end-1:end), 'am') % for '_am'
                    gasfolder_short = gasfolder{i}(1:end-3);
                elseif length(gasfolder{i}) > 3 && strcmp(gasfolder{i}(end-1:end), 'pm') % for '_pm'
                    gasfolder_short = gasfolder{i}(1:end-3);
                else
                    gasfolder_short = gasfolder{i};
                end
                % create and add the data to the file
                ncfilename_j = sprintf('SPARC_DI_T2Mz_%s_%i_ACEFTS_v3.6_i01.nc', gasfolder{i}, data_years_unique(j));
                outdir = fullfile(ncdir, gasfolder_short);
                if exist(outdir,'dir') ~= 7
                    mkdir(outdir);
                end
                ncfilename_j = fullfile(outdir, ncfilename_j);
                make_ace_empty_netcdf_SDI(ncfilename_j,gasfolder_short); % make the empty .nc file (no data)
                fprintf('writing the %i data...\n', data_years_unique(j))
                ncwrite(ncfilename_j, gasfolder_short, vmr_out);
                ncwrite(ncfilename_j, strcat(gasfolder_short,'_STD'), vmr_std_out);
                ncwrite(ncfilename_j, strcat(gasfolder_short,'_NR'), vmr_obs_out);
                ncwrite(ncfilename_j, 'LST_MEAN', lst_mean);
                ncwrite(ncfilename_j, 'LST_MAX', lst_max);
                ncwrite(ncfilename_j, 'LST_MIN', lst_min);
                ncwrite(ncfilename_j, 'AVE_DOM', ave_dom);
                ncwrite(ncfilename_j, 'AVE_LAT', ave_lat);
                ncwrite(ncfilename_j, 'lat', lat);
                ncwrite(ncfilename_j, 'plev', plev);
                ncwrite(ncfilename_j, 'time', time);
                disp('Done')
            end
        end
    end 
else
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',ncdir)
end
disp('All done :)')
%
end


function [ ] = make_ace_relationship_season( tanstruct_x_in, tanstruct_y_in, out_directory)
%A function to create zonally averaged climatologies of ACE measurements,
%by each unique calendar month. 'make_ace_climatology.m' is called here.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. The GLC data must also be added to
%           the tanstruct so that it has the latitude information.
%
%           out_directory: STRING - the path to the directory in which you
%           would like the output to be saved.
%
% *OUTPUT*
%           .mat files of STRUCTURES will be written to 'out_directory'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 09/18

tic
%% Things that may be changed often
if nargin < 2
    home_linux = '/home/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
    home_mac = '/Users/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
    home_windows = 'C:\Users\ryann\Dropbox\climatology\nryan\'; %#ok<NASGU>
    reldirectory = strcat(home_mac,'climdata/');
else
    reldirectory = strcat(out_directory,'/');
end

%the name of the output files
% savename_pre = 'CMAM_CLIM_MLSsample_lat_';
savename_pre = 'ACEFTS_REL_v3_lat_';
% cells with the names of the month
monthnames = {'January', 'February', 'March', 'April', 'May', 'June', ...
    'July', 'August', 'September', 'October', 'November', 'December'};

%% Define some things
gasx = tanstruct_x_in;
gasy = tanstruct_y_in;
relstruct = []; %#ok<NASGU>
gasname_out = strcat(gasx.gas,'-',gasy.gas);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lat1 = 0;
lat2 = 30;
gasx = subset_ace_by_lat(gasx,lat1,lat2);
gasy = subset_ace_by_lat(gasy,lat1,lat2);
gasname_out = strcat(gasname_out,'-','0030N');

%% Only continue of the output directory exists
if isdir(reldirectory)
    %% Create a folder to store the climatology data if one does not exist
    reldirectory_gas = strcat(reldirectory,gasname_out,'/','serial_month','/')
    if exist(reldirectory_gas,'dir') ~= 7
        mkdir(reldirectory_gas);
    end
    %% match the data
    [gasx,gasy] = match_ace_data(gasx,gasy);
    if isempty(gasx.occultation)
        error('there is no overlap between the input datasets')
    end
    % get the unique years for the data
    vdates = datevec(mjd2datenum(gasx.date_mjd));
    years_unique = unique(vdates(:,1));
    
    %% For each year loop through the months and create climatologies for each.
    for j = 1:length(years_unique)
        %subset the data by year
        warning off % suppress warnings about reducing the data to zero here. There is output below if this is the case
        gasx_yearj = subset_ace_by_year(gasx,years_unique(j));
        gasy_yearj = subset_ace_by_year(gasy,years_unique(j));
        warning on
        for i = 1:12
            %subset the ace data by month
            warning off % suppress warnings about reducing the data to zero here. There is output below if this is the case
            gasx_yearj_monthi = subset_ace_by_month(gasx_yearj,i);
            gasy_yearj_monthi = subset_ace_by_month(gasy_yearj,i);
            warning on
            if isempty(gasx_yearj_monthi.date_mjd) == 0
                fprintf('\nPreparing relationship for %s %d', monthnames{i}, years_unique(j))
                relstruct_monthi = make_ace_relationship(gasx_yearj_monthi, gasy_yearj_monthi);
                relstruct_monthi.climatology_type = 'serial_month';
                relstruct_monthi.time = [years_unique(j),i];
                
                % save the file as a matlab structure for now. In the chosen directory
                if  ~isempty(relstruct_monthi.date_mjd)
                    relstruct = relstruct_monthi; %#ok<NASGU>
                    savename_post = sprintf('_%d_%02d',years_unique(j),i);
                    savedest = strcat(reldirectory_gas,savename_pre,gasname_out,savename_post);
                    fprintf('\nSaving %s gas relation to %s\n', monthnames{i}, savedest);
                    save(savedest,'relstruct');
                else
                    fprintf('There are no data for %s %d\n', monthnames{i}, years_unique(j));
%                     j
                    years_unique(j)
%                     i
                end
            else
                fprintf('There are no ACE data for %s %d\n', monthnames{i}, years_unique(j));
            end
        end
    end
else
    fprintf('\nThe directory, %s, does not exist. Please create it first\n', reldirectory)
end
fprintf('\nDone :)\n')
toc
%
end
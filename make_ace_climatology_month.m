function [ ] = make_ace_climatology_month( tanstruct, out_directory)
%A function to create zonally averaged climatologies of ACE measurements,
%by calendar month. 'make_ace_climatology.m' is called here.

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
%   NJR - 02/18

tic
%% Things that may be changed often
if nargin < 2
    home_linux = '/home/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
    home_mac = '/Users/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
    home_windows = 'C:\Users\ryann\Dropbox\climatology\nryan\'; %#ok<NASGU>
    climdirectory = strcat(home_windows,'climdata/');
else
    climdirectory = strcat(out_directory,'/');
end

%the name of the output files
savename_pre = 'ACEFTS_CLIM_v3_lat_';
% cells with the names of the month
monthnames = {'January', 'February', 'March', 'April', 'May', 'June', ...
    'July', 'August', 'September', 'October', 'November', 'December'};

%% Define some things
gas = tanstruct;
climstruct = []; %#ok<NASGU>

%% Create a folder to store the climatology data if one does not exist
climdirectory_gas = strcat(climdirectory,gas.gas,'/');
if exist(climdirectory_gas,'dir') ~= 7
    mkdir(climdirectory_gas);
end

%% Only continue of the output directory exists
if isdir(climdirectory)  
    %% Loop through the months and create climatologies for each.
    for i = 1:12
        %subset the ace data by month
        warning off % supress warnings about reducing the data to zero here. There is output below if this is the case
        gas_monthi = subset_ace_by_month(gas,i);
        warning on
        fprintf('\nPreparing climatology for %s', monthnames{i})
        climstruct_monthi = make_ace_climatology(gas_monthi);
        climstruct_monthi.climatology_type = 'calendar_month';
        climstruct_monthi.time = i;
        
        % save the file as a matlab structure for now. In the chosen directory
        if  ~isempty(climstruct_monthi.start_date)
            climstruct = climstruct_monthi; %#ok<NASGU>
            savename_post = sprintf('_%02d',i);
            savedest = strcat(climdirectory_gas,savename_pre,gas.gas,savename_post);
            fprintf('\nSaving %s climatology to %s\n', monthnames{i}, savedest);
            save(savedest,'climstruct');
        else
            fprintf('There are no climatology data for %s\n', monthnames{i});
        end
    end
else
    fprintf('\nThe directory, %s, does not exist. Please create it first\n', climdirectory)
end
fprintf('\nDone :)\n')
toc
%
end
function [ ] = make_ace_climatology_3month_by_eql( tanstruct, out_directory)
%A function to create zonally averaged climatologies of ACE measurements
%by 3-monthly periods, using equivalent latitude.
%'make_ace_climatology_by_eql.m' is called here. 

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. The DMP data must also be added to
%           the tanstruct so that it has the equivalent latitude
%           information.
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
savename_pre = 'ACEFTS_CLIM_v3_eql_';
% cells with the names of the month
monthnames = {'DJF', 'MAM', 'JJA', 'SON'};

%% Define some things
gas = tanstruct;
climstruct = []; %#ok<NASGU>

%% Only continue of the output directory exists
if isdir(climdirectory)
    %% Create a folder to store the climatology data if one does not exist
    climdirectory_gas = strcat(climdirectory,gas.gas,'/');
    if exist(climdirectory_gas,'dir') ~= 7
        mkdir(climdirectory_gas);
    end
    %% Loop through the months and create climatologies for each.
    for i = 1:4
        %subset the ace data by month
        warning off % supress warnings about reducing the data to zero here. There is output below if this is the case
        gas_3monthi = subset_ace_by_3month(gas,monthnames{i});
        warning on
        fprintf('\nPreparing climatology for %s...', monthnames{i})
        climstruct_3monthi = make_ace_climatology_by_eql(gas_3monthi);
        climstruct_3monthi.climatology_type = 'season';
        climstruct_3monthi.time = i;
        
        % save the file as a matlab structure for now. In the chosen directory
        if  ~isempty(climstruct_3monthi.start_date)
            climstruct = climstruct_3monthi; %#ok<NASGU>
            savename_post = sprintf('_%s',monthnames{i});
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
function [ ] = make_ace_climatology_multiple_by_eql( clim_type, varargin )
%A function to read theACE v3.5/6 .mat data and calculate a climatology,
%using equivalent latitude, with the relevant information. The directory
%containing the .nc data, as well as the directory to which you want to
%write the .mat data, should be defined below. The directory containing
%the ACE .mat files and the output directory for the climatologies are
%defined by the user below.

% *INPUT*
%           clim_type: STRING - the type of climatology that you would like
%           to make. Currently, the options are
%           month        a climatology averaged by calendar month. Produces
%                        12 files. Calls 'make_ace_climatology_month.m'.
%           3month       a climatology averaged by 3-month periods of DJF,
%                        MAM, JJA, and SON. Produces 4 files. Calls
%                        'make_ace_climatology_3month.m'.
%           serialmonth  a climatology of each unique calendar month in the
%                        data. Produces as many files as there are unique
%                        calendar months in the data. Calls
%                        'make_ace_climatology_serialmonth.m'.
%
%           varargin: STRING - the name of the gas for which you want to
%                calculate a climatology. The directory containing the .nc
%                data, as well as the directory to which you want to write
%                the .mat data, should be defined below. The ACE GLC data
%                is also read here and used for the climatology.
%
% *OUTPUT*
%                .mat files will be written to the output folder that is
%                specified below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 02/18

%% Define some things
% the naming convention of ACE .mat files is ACE_v3p6_gas.mat, where
% 'gas' is the gas species. e.g., ACE_v3p6_O3.mat

%%USER DEFINED
home_linux = '/home/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\Dropbox\climatology\nryan\'; %#ok<NASGU>
home_deluge = '/net/deluge/pb_1/users/nryan/';

% matdirectory = '/Volumes/Seagate Backup Plus Drive/ACE/matdata/';
% matdirectory = 'F:\ACE\matdata\';
% matdirectory = strcat(home_mac,'matdata/');
matdirectory = strcat(home_deluge,'ACE/','matdata/');
if ~isdir(matdirectory)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',matdirectory)
    error('The directory containing the .nc data couldn''t be found')
end
% climdirectory = '/Users/niall/Dropbox/climatology/nryan/matclim/'; % edit this to your output directory
% climdirectory = '/Volumes/Seagate Backup Plus Drive/ACE/matclim/';
% climdirectory = 'F:\ACE\climdata\';
% climdirectory = strcat(home_mac,'climdata/');
climdirectory = strcat(home_deluge,'ACE/','climdata/');

%%STANDARD
filein_pre = 'ACE_v3p6_';
filein_post = '.mat';

%% Read the DMP file
filein = strcat(matdirectory,filein_pre,'DMPv2p0',filein_post);
temp_dir = dir(filein); % to check whether the GLC file exists
if(~isempty(temp_dir)) % if the file exists
    dmp = load(filein); dmp = dmp.dmpstruct;
else % if the file doesn't exist
    error('Error: I can''t find %s\n',filein);
end

%% Read and save the data
if isdir(climdirectory)
    
    %% Get the names of the input gases
    gases = varargin; % cell with gas names
    lin = length(gases); % number of gases to read
    climtype = clim_type;
    
    %% Read the files
    if lin == 1
        gasin = gases{1};
        if ~strcmp(gasin,'all') %%% a case of reading one .mat file
            filein = strcat(matdirectory,filein_pre,gasin,filein_post);
            temp_dir = dir(filein); % to check whether the file exists
            if(~isempty(temp_dir)) % if the file exists
                fprintf('\nPROCESSING %s\n',gasin)
                gasi = load(filein); gasi = gasi.tanstruct;
                gasidmp = merge_ace_dmp(gasi,dmp);
                switch climtype
                    case 'month'
                        make_ace_climatology_month_by_eql(gasidmp,climdirectory);
                    case '3month'
                        make_ace_climatology_3month_by_eql(gasidmp,climdirectory);
                    case 'serialmonth'
                        make_ace_climatology_serialmonth_by_eql(gasidmp,climdirectory);
                end
            else
                fprintf('Error: I can''t find %s\n',filein); % if the file doesn't exist
                fprintf('moving on...\n')
            end
        elseif strcmp(gasin,'all') %%% a case of wanting to read all of the .nc files in the ncdirectory
            tempdir = dir(strcat(matdirectory,'ACE_v3p6_*.mat')); % gets the info about the ACE .mat files in the ncdirectory
            fileall = {tempdir.name}; % a cell of the filenames
            %             fileall{1}(10:end-4)
            if ~isempty(fileall)
                for i = 1:length(fileall)
                    filein = strcat(matdirectory,fileall{i}); % get the path to the ith file in the list
                    if ~strcmp(fileall{i}(10:end-4), 'DMPv2p0') && ~strcmp(fileall{i}(10:end-4), 'GLC')
                        fprintf('\nPROCESSING %s\n',fileall{i}(10:end-4))
                        gasi = load(filein); gasi = gasi.tanstruct;
                        gasidmp = merge_ace_dmp(gasi,dmp);
                        switch climtype
                            case 'month'
                                make_ace_climatology_month_by_eql(gasidmp,climdirectory);
                            case '3month'
                                make_ace_climatology_3month_by_eql(gasidmp,climdirectory);
                            case 'serialmonth'
                                make_ace_climatology_serialmonth_by_eql(gasidmp,climdirectory);
                        end
                        clear gasidmp
                    end
                end
            else
                fprintf('I couldn''t find any ACE .nc data in %s.\nIs the path correct?\n',ncdirectory)
            end
        end
    else % if there are multiple gas names as the input
        for i = 1:lin
            gasin = gases{i};
            filein = strcat(matdirectory,filein_pre,gasin,filein_post);
            temp_dir = dir(filein); % to check whether the file exists
            if(~isempty(temp_dir)) % if the file exists
                fprintf('\nPROCESSING %s\n',gasin)
                gasi = load(filein); gasi = gasi.tanstruct;
                gasidmp = merge_ace_dmp(gasi,dmp);
                switch climtype
                    case 'month'
                        make_ace_climatology_month_by_eql(gasidmp,climdirectory);
                    case '3month'
                        make_ace_climatology_3month_by_eql(gasidmp,climdirectory);
                    case 'serialmonth'
                        make_ace_climatology_serialmonth_by_eql(gasidmp,climdirectory);
                end
                clear gasiglc
            else
                fprintf('Error: I can''t find %s\n',filein); % if the file doesn't exist
                fprintf('moving on...\n')
            end
        end
    end
    
else
    fprintf('\nThe directory, %s, does not exist. Please create it first\n', climdirectory)
end
%
end


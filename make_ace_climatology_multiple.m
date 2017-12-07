function [ ] = make_ace_climatology_multiple( varargin )
%A function to read theACE v3.5/6 .mat data and calculate a climatology
%with the relevant information. The directory containing the .nc data, as
%well as the directory to which you want to write the .mat data, should be
%defined below. The ACE GLC data is also read here and saved as a .mat
%structure. The direcetory containing the ACE data .mat files and the
%output directory for the climatologies is defined by the user below

% *INPUT*
%           gas: STRING - the name of the gas for which you want to
%                calculate a climatology. The directory containing the .nc
%                data, as well as the directory to which you want to write
%                the .mat data, should be defined below. The ACE GLC data
%                is also read here and used for the climatology.
%
% *OUTPUT*
%                .mat files will be written to the outut folder that is
%                specified below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 11/17

%% Define some things
% the naming convention of ACE .mat files is ACE_v3p6_gas.mat, where
% 'gas' is the gas species. e.g., ACE_v3p6_O3.mat

%%USER DEFINED
% matdirectory = '/Users/niall/Dropbox/climatology/nryan/matdata/'; % edit this to your directory that contains the ACE .mat data
matdirectory = '/Volumes/Seagate Backup Plus Drive/ACE/matdata/';
if ~isdir(matdirectory)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',matdirectory)
    error('The directory containing the .nc data couldn''t be found')
end
% climdirectory = '/Users/niall/Dropbox/climatology/nryan/matclim/'; % edit this to your output directory
climdirectory = '/Volumes/Seagate Backup Plus Drive/ACE/matclim/';

%%STANDARD
filein_pre = 'ACE_v3p6_';
filein_post = '.mat';

%% Read the GLC file
filein = strcat(matdirectory,filein_pre,'GLC',filein_post);
temp_dir = dir(filein); % to check whether the GLC file exists
if(~isempty(temp_dir)) % if the file exists
    glc = load(filein); glc = glc.glcstruct;
else % if the file doesn't exist
    error('Error: I can''t find %s\n',filein);
end

%% Read and save the data
if isdir(climdirectory)
    
    %% Get the names of the input gases
    gases = varargin; % cell with gas names
    lin = length(gases); % number of gases to read
    
    %% Read the files
    if lin == 1
        gasin = gases{1};
        if ~strcmp(gasin,'all') %%% a case of reading one .mat file
            filein = strcat(matdirectory,filein_pre,gasin,filein_post);
            temp_dir = dir(filein); % to check whether the file exists
            if(~isempty(temp_dir)) % if the file exists
                fprintf('\nPROCESSING %s\n',gasin)
                gasi = load(filein); gasi = gasi.tanstruct;
                gasiglc = merge_ace_glc(gasi,glc);
                make_ace_climatology(gasiglc,climdirectory);
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
                    if ~strcmp(fileall{i}(10:end-4), 'GLC')
                        fprintf('\nPROCESSING %s\n',fileall{i}(10:end-4))
                        gasi = load(filein); gasi = gasi.tanstruct;
                        gasiglc = merge_ace_glc(gasi,glc);
                        make_ace_climatology(gasiglc,climdirectory);
                        clear gasiglc
                    end
                end
            else
                fprintf('I couldn''t find any ACE .nc data in %s.\nIs the path correct?\n',ncdirectory)
            end
        end
    else % if there are multiple gas names input
        for i = 1:lin
            gasin = gases{i};
            filein = strcat(matdirectory,filein_pre,gasin,filein_post);
            temp_dir = dir(filein); % to check whether the file exists
            if(~isempty(temp_dir)) % if the file exists
                fprintf('\nPROCESSING %s\n',gasin)
                gasi = load(filein); gasi = gasi.tanstruct;
                gasiglc = merge_ace_glc(gasi,glc);
                make_ace_climatology(gasiglc,climdirectory);
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


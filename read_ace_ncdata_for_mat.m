function [ ] = read_ace_ncdata_for_mat( varargin )
%A function to read theACE v3.5/6 .nc data and create .mat files with
%relevant information. The directory containing the .nc data, as well as
%the directory to which you want to write the .mat data, should be defined
%below. The ACE GLC data is also read here and saved as a .mat structure.

% *INPUT*
%           gas: STRING - the name of the gas for which you want to read
%                the data. Use 'GLC' to read the GLC data.
%                The input may be a single gas (e.g., 'O3'), or multiple
%                gases (e.g., 'O3', 'ClO').
%                To read all .nc files in a directory, the input is 'all'.
%
% *OUTPUT*
%                .mat files will be written to the outut folder that is
%                specified below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 10/17

%% Define some things
% the naming convention of ACE .nc files is ACEFTS_L2_v3p6_gas.nc, where
% 'gas' is the gas species. e.g., ACEFTS_L2_v3p6_O3.nc

%%USER DEFINED
% ncdirectory = '/Users/niall/Dropbox/climatology/nryan/ncdata/'; % edit this to your directory that contains the ACE netcdf data
% ncdirectory = '/Volumes/Seagate Backup Plus Drive/ACE/ncdata/v3p6/';
% ncdirectory = 'F:\ACE\climdata\';
ncdirectory = '/net/deluge/pb_1/users/nryan/ACE/ncdata/v3p6/';
if ~isdir(ncdirectory)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',ncdirectory)
    error('The directory containing the .nc data couldn''t be found')
end
% matdirectory = '/Users/niall/Dropbox/climatology/nryan/matdata/'; % edit this to your output directory
% matdirectory = '/Volumes/Seagate Backup Plus Drive/ACE/matdata/';
% matdirectory = 'F:\ACE\matdata\';
matdirectory = '/net/deluge/pb_1/users/nryan/ACE/matdata/';

%%STANDARD
filein_pre = 'ACEFTS_L2_v3p6_';
filein_post = '.nc';
savename_pre = 'ACE_v3p6_';

%% Read and save the data
if isdir(matdirectory)
    
    %% Get the names of the input gases
    gases = varargin; % array of cells with gas names
    lin = length(gases); % number of gases to read
    
    %% Read the files
    if lin == 1
        gasin = gases{1};
        if ~strcmp(gasin,'all') %%% a case of reading one .nc file
            filein = strcat(ncdirectory,filein_pre,gasin,filein_post);
            temp_dir = dir(filein); % to check whether the file exists
            if(~isempty(temp_dir)) % if the file exists
                if strcmp(gasin, 'GLC')
                    [glcstruct, gasouti] = read_ace_ncdata_glc(filein); %#ok<ASGLU> % read the data and create a matlab structure
                    savedest = strcat(matdirectory,savename_pre,gasouti);
                    fprintf('saving %s data to %s\n', gasouti, savedest);
                    save(savedest,'glcstruct');
                    fprintf('done\n')
                    clear glcstruct
                else
                    [tanstruct, gasouti] = read_ace_ncdata(filein); %#ok<ASGLU> % read the data and create a matlab structure
                    savedest = strcat(matdirectory,savename_pre,gasouti);
                    fprintf('saving %s data to %s\n', gasouti, savedest);
                    save(savedest,'tanstruct');
                    fprintf('done\n')
                end
                savedest = strcat(matdirectory,savename_pre,gasouti);
                fprintf('saving %s data to %s\n', gasouti, savedest);
                save(savedest,'tanstruct');
                fprintf('done\n')
            else
                fprintf('Error: I can''t find %s\n',filein); % if the file doesn't exist
                fprintf('moving on...\n')
            end
        elseif strcmp(gasin,'all') %%% a case of wanting to read all of the .nc files in the ncdirectory
            tempdir = dir(strcat(ncdirectory,'ACEFTS_L2_v3p6*.nc')); % gets the info about the ACE .nc files in the ncdirectory
            fileall = {tempdir.name}; % a cell of the filenames
            %fileall{1}(16:end-3)
            if ~isempty(fileall)
                for i = 1:length(fileall)
                    filein = strcat(ncdirectory,fileall{i}); % get the path to the ith file in the list
                    if strcmp(fileall{i}(16:end-3), 'GLC')
                        [glcstruct, gasouti] = read_ace_ncdata_glc(filein); %#ok<ASGLU> % read the data and create a matlab structure
                        savedest = strcat(matdirectory,savename_pre,gasouti);
                        fprintf('saving %s data to %s\n', gasouti, savedest);
                        save(savedest,'glcstruct');
                        fprintf('done\n')
                        clear glcstruct
                    else
                        [tanstruct, gasouti] = read_ace_ncdata(filein); %#ok<ASGLU> % read the data and create a matlab structure
                        savedest = strcat(matdirectory,savename_pre,gasouti);
                        fprintf('saving %s data to %s\n', gasouti, savedest);
                        save(savedest,'tanstruct');
                        fprintf('done\n')
                    end
                end
            else
                fprintf('I couldn''t find any ACE .nc data in %s.\nIs the path correct?\n',ncdirectory)
            end
        end
    else % if there are multiple gas names input
        for i = 1:lin
            gasin = gases{i};
            filein = strcat(ncdirectory,filein_pre,gasin,filein_post);
            temp_dir = dir(filein); % to check whether the file exists
            if(~isempty(temp_dir)) % if the file exists
                if strcmp(gasin, 'GLC')
                    [glcstruct, gasouti] = read_ace_ncdata_glc(filein); %#ok<ASGLU> % read the data and create a matlab structure
                    savedest = strcat(matdirectory,savename_pre,gasouti);
                    fprintf('saving %s data to %s\n', gasouti, savedest);
                    save(savedest,'glcstruct');
                    fprintf('done\n')
                    clear glcstruct
                else
                    [tanstruct, gasouti] = read_ace_ncdata(filein); %#ok<ASGLU> % read the data and create a matlab structure
                    savedest = strcat(matdirectory,savename_pre,gasouti);
                    fprintf('saving %s data to %s\n', gasouti, savedest);
                    save(savedest,'tanstruct');
                    fprintf('done\n')
                end
                savedest = strcat(matdirectory,savename_pre,gasouti);
                fprintf('saving %s data to %s\n', gasouti, savedest);
                save(savedest,'tanstruct');
                fprintf('done\n')
            else
                fprintf('Error: I can''t find %s\n',filein); % if the file doesn't exist
                fprintf('moving on...\n')
            end
        end
    end
    
else
    fprintf('\nThe directory, %s, does not exist. Please create it first\n', matdirectory)
end
%
end


function [ ] = write_ace_scaled_apriori( varargin )
%A function to read the ACE v3.5/6 .mat data and alter the quality flags in
%the structure so that the flags indicating scaled a priori information in
%measurements are replaced by flags that indicate ordinary mesurement data.
%A new file is created that contains the information. The name of the
%species in the output file (and output filename) is changed to 'GAS_sap',
%with GAS indicating the original name of the species.
%
% *INPUT*
%
%           varargin: STRING - the name of the gas for which you want to
%                make files that treat scaled apriori data as measurements.
%                The directory containing the .mat data, which is the
%                directory to which you want to write the output files,
%                should be defined below. 
%
% *OUTPUT*
%                .mat files will be written to the output folder that is
%                specified below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 03/18

%% Define some things
% the naming convention of ACE .mat files is ACE_v3p6_gas.mat, where
% 'gas' is the gas species. e.g., ACE_v3p6_O3.mat

%%USER DEFINED
home_linux = '/home/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\Dropbox\climatology\nryan\'; %#ok<NASGU>
home_deluge = '/net/deluge/pb_1/users/nryan/'; %#ok<NASGU>

matdirectory = 'C:\Users\ryann\ACE\matdata\';
% matdirectory = '/Volumes/Seagate Backup Plus Drive/ACE/matdata/';
% matdirectory = 'F:\ACE\matdata\';
% matdirectory = strcat(home_mac,'matdata/');
% matdirectory = strcat(home_deluge,'ACE/','matdata/');
if ~isdir(matdirectory)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',matdirectory)
    error('The directory containing the .mat data couldn''t be found')
end

%%STANDARD
filein_pre = 'ACE_v3p6_';
filein_post = '.mat';
fileout_pre = filein_pre;

%% Get the names of the input gases
    gases = varargin; % cell with gas names
    lin = length(gases); % number of gases to read
    
    if lin == 1
        gasin = gases{1};
        if ~strcmp(gasin,'all') %%% a case of reading one .mat file
            filein = strcat(matdirectory,filein_pre,gasin,filein_post);
            temp_dir = dir(filein); % to check whether the file exists
            if(~isempty(temp_dir)) % if the file exists
                fprintf('\nPROCESSING %s\n',gasin)
                gasi = load(filein); gasi = gasi.tanstruct;
                gasi_sap = include_ace_scaled_apriori(gasi);
                % save the data
                savedest = strcat(matdirectory,fileout_pre,gasin,'_sap');
                fprintf('saving %s_sap data to %s\n', gasin, savedest);
                tanstruct = gasi_sap;
%                 tanstruct.gas = strcat(gasin,'_sap'); % this is already
%                 done in the above function
                save(savedest,'tanstruct');
                fprintf('done\n')
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
                    if ~strcmp(fileall{i}(10:end-4), 'GLC') && ~strcmp(fileall{i}(10:end-4), 'DMPv2p0')
                        fprintf('\nPROCESSING %s\n',fileall{i}(10:end-4))
                        gasi = load(filein); gasi = gasi.tanstruct;
                        gasi_sap = include_ace_scaled_apriori(gasi);
                        % save the am data
                        savedest = strcat(matdirectory,fileout_pre,gasin,'_sap');
                        fprintf('saving %s_sap data to %s\n', gasin, savedest);
                        tanstruct = gasi_sap;
%                         tanstruct.gas = strcat(gasin,'_sap'); % this is already
%                           done in the above function
                        save(savedest,'tanstruct');
                        fprintf('done\n')
                        clear tanstruct gasi
                    end
                end
            else
                fprintf('I couldn''t find any ACE v3p6 .mat data in %s.\nIs the path correct?\n',matdirectory)
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
                gasi_sap = include_ace_scaled_apriori(gasi);
                % save the data
                savedest = strcat(matdirectory,fileout_pre,gasin,'_sap');
                fprintf('saving %s_sap data to %s\n', gasin, savedest);
                tanstruct = gasi_sap; 
%                 tanstruct.gas = strcat(gasin,'_sap'); % this is already
%                 done in the above function
                save(savedest,'tanstruct');
                fprintf('done\n')
                clear tanstruct gasi
            else
                fprintf('Error: I can''t find %s\n',filein); % if the file doesn't exist
                fprintf('moving on...\n')
            end
        end
    end


end


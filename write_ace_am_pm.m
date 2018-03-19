function [ ] = write_ace_am_pm( varargin )
%A function to read th eACE v3.5/6 .mat data and split the data by am and
%pm measurements. am and pm are determined by local solar time (LST).

% *INPUT*
%
%           varargin: STRING - the name of the gas for which you want to
%                make am/pm files. The directory containing the .mat
%                data, which is the directory to which you want to write
%                the output files, should be defined below. The ACE GLC
%                data is also read here and used for determining the LST.
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

matdirectory = '/Volumes/Seagate Backup Plus Drive/ACE/matdata/';
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

%% Read the GLC file
filein = strcat(matdirectory,filein_pre,'GLC',filein_post);
temp_dir = dir(filein); % to check whether the GLC file exists
if(~isempty(temp_dir)) % if the file exists
    glc = load(filein); glc = glc.glcstruct;
else % if the file doesn't exist
    error('Error: I can''t find %s\n',filein);
end

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
                gasiglc = merge_ace_glc(gasi,glc);
                [gasi_am, gasi_pm] = subset_ace_by_lst(gasiglc); % no argument uses noon as the split time
                % save the am data
                savedest_am = strcat(matdirectory,fileout_pre,gasin,'_am');
                fprintf('saving %s_am data to %s\n', gasin, savedest_am);
                tanstruct = gasi_am;
                tanstruct.gas = strcat(gasin,'_am');
                save(savedest_am,'tanstruct');
                %save the pm data
                savedest_pm = strcat(matdirectory,fileout_pre,gasin,'_pm');
                fprintf('saving %s_pm data to %s\n', gasin, savedest_pm);
                tanstruct = gasi_pm;
                tanstruct.gas = strcat(gasin,'_pm');
                save(savedest_pm,'tanstruct');
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
                        gasiglc = merge_ace_glc(gasi,glc);
                        [gasi_am, gasi_pm] = subset_ace_by_lst(gasiglc); % no argument uses noon as the split time
                        % save the am data
                        savedest_am = strcat(matdirectory,fileout_pre,gasin,'_am');
                        fprintf('saving %s_am data to %s\n', gasin, savedest_am);
                        tanstruct = gasi_am;
                        tanstruct.gas = strcat(gasin,'_am');
                        save(savedest_am,'tanstruct');
                        %save the pm data
                        savedest_pm = strcat(matdirectory,fileout_pre,gasin,'_pm');
                        fprintf('saving %s_pm data to %s\n', gasin, savedest_pm);
                        tanstruct = gasi_pm;
                        tanstruct.gas = strcat(gasin,'_pm');
                        save(savedest_pm,'tanstruct');
                        fprintf('done\n')
                        clear tanstruct gasiglc
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
                gasiglc = merge_ace_glc(gasi,glc);
                [gasi_am, gasi_pm] = subset_ace_by_lst(gasiglc); % no argument uses noon as the split time
                % save the am data
                savedest_am = strcat(matdirectory,fileout_pre,gasin,'_am');
                fprintf('saving %s_am data to %s\n', gasin, savedest_am);
                tanstruct = gasi_am; 
                tanstruct.gas = strcat(gasin,'_am');
                save(savedest_am,'tanstruct');
                %save the pm data
                savedest_pm = strcat(matdirectory,fileout_pre,gasin,'_pm');
                fprintf('saving %s_pm data to %s\n', gasin, savedest_pm);
                tanstruct = gasi_pm;
                tanstruct.gas = strcat(gasin,'_pm');
                save(savedest_pm,'tanstruct');
                fprintf('done\n')
                clear tanstruct gasiglc
            else
                fprintf('Error: I can''t find %s\n',filein); % if the file doesn't exist
                fprintf('moving on...\n')
            end
        end
    end


end


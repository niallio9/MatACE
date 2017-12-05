function [ ] = read_ace_dmp_for_mat( path_to_dmp_data )
%A function to read the ACE v3.5/6 DMP v2.0 files and make a .mat
%structure. This function assumes that all of the DMP files are in the same
%directory. The script 'move_files_up_2levels.m' can be used to move files
%from the original file structure to a single folder. The output folder for
%the .mat data should be specified by the user below.

% *INPUT*
%           path_to_dmp_data: STRING - the path to the directory that holds
%           the dmp files.The script 'move_files_up_2levels.m' can be
%           used to move files from the original file structure to a
%           single folder. 
%
% *OUTPUT*
%           .mat files will be written to the outut folder that is
%           specified below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 10/2017

%% Define some things
dmpdirectory = path_to_dmp_data;
if ~isdir(path_to_dmp_data)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',path_to_dmp_data)
    error('The directory containing the dmp files couldn''t be found')
end

%%USER DEFINED
matdirectory = '/Users/niall/Dropbox/climatology/nryan/dmptest10'; % edit this to your output directory

%%STANDARD
savename = 'ACE_v3p5_6_DMP2p0';

%%
if isdir(matdirectory) 
    %% get a list of the .asc files in the directory
    tempdir = dir(strcat(dmpdirectory,'/','*.asc'));
    fileall = {tempdir.name};
    lfiles = length(fileall);
    
    
    %% Read out the data from each of the files and put them into a structure
    % There seems to be 76 levels of dmp data for now.
    rowdmp = 76;
    coldmp = lfiles;
    out.source_file = cell(1,coldmp); %string
    out.occultation = nan(1,coldmp);
    out.sr1ss0 = nan(1,coldmp);
    out.date_mjd = nan(1,coldmp);
    out.lon_tangent = nan(1,coldmp);
    out.lat_tangent  = nan(1,coldmp);
    out.version = cell(1,coldmp); % string
    out.altitude_km = nan(rowdmp,1); % all the altitude grids are the same
    out.T = nan(rowdmp,coldmp);
    out.pressure_hPa = nan(rowdmp,coldmp);
    out.lon = nan(rowdmp,coldmp);
    out.lat = nan(rowdmp,coldmp);
    out.Theta = nan(rowdmp,coldmp);
    out.spv = nan(rowdmp,coldmp);
    out.eql = nan(rowdmp,coldmp);
    
    fprintf('\nReading files...\n')
    for i = 1:lfiles
        filei = strcat(dmpdirectory,'/',fileall{i});
        fprintf('%s\n',filei);
        [ dmpi ] = read_ace_dmpv2(filei);
        out.source_file{i} = dmpi.source_file;
        out.occultation(i) = dmpi.occultation;
        out.sr1ss0(i) = dmpi.sr1ss0;
        out.date_mjd(i) = dmpi.date_mjd;
        out.lon_tangent(i) = dmpi.lon_tangent;
        out.lat_tangent(i) = dmpi.lat_tangent;
        out.version{i} = dmpi.version;
        if i == 1
            out.altitude_km(:,i) = dmpi.altitude_km; % all the altitude grids are the same
        end
        out.T(:,i) = dmpi.T;
        out.pressure_hPa(:,i) = dmpi.pressure_hPa;
        out.lon(:,i) = dmpi.lon;
        out.lat(:,i) = dmpi.lat;
        out.Theta(:,i) = dmpi.Theta;
        out.spv(:,i) = dmpi.spv;
        out.eql(:,i) = dmpi.eql;
    end
    dmpstruct = out; %#ok<NASGU>
    
    savedest = strcat(matdirectory,'/',savename);
    fprintf('\nsaving DMP data to %s\n', savedest);
    save(savedest,'dmpstruct');
    fprintf('done\n')
else
    fprintf('\nThe directory, %s, does not exist. Please create it first\n', matdirectory)
end
end

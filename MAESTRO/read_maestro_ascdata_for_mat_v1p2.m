function [ tanstruct ] = read_maestro_ascdata_for_mat( path_to_maestro_data )
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
datadirectory = path_to_maestro_data;
if ~isdir(path_to_maestro_data)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',path_to_maestro_data)
    error('The directory containing the dmp files couldn''t be found')
end

%%USER DEFINED
% matdirectory = '/Users/niall/Dropbox/climatology/nryan/dmptest10'; % edit this to your output directory
matdirectory = 'C:\Users\ryann\ACE\MAESTRO\matdata'; % edit this to your output directory
% matdirectory = 'F:\ACE\dmp3.5';

%%STANDARD
savename = 'MAESTRO_v1p2_O3u';

%%
if isdir(matdirectory) 
    
    %% get a list of the folders in the directory
    tempdir = dir(datadirectory);
    datafolders = {tempdir.name}; % asuming that all that is in the directory are the folders with the data
    datafolders = datafolders(3:end); % get rid of ',' and'..' 
    
    %% Read out the data from each of the files and put them into a structure
    % There seems to be 311 levels of O3 data for now.
    rownum = 201;
    colnum = 2e5; % i think there's about 1.2e5 measurements from maestro at the moment, so this is an overshoot
    out.source_file = cell(1,colnum); %string
    out.version = cell(1,colnum); % string
%     out.source_file = cell(1,1);
%     out.version = cell(1,1);
    out.occultation = nan(1,colnum);
    out.sr1ss0 = nan(1,colnum);
    out.beta_angle = nan(1,colnum);
    out.date_mjd = nan(1,colnum);
    out.gas ='O3v'; % this is fexed for now as I'm only doing O3 at the moment
    out.altitude_km = nan(rownum,colnum ); % all the altitude grids are the same
    out.vmr = nan(rownum,colnum);
    out.vmr_error = nan(rownum,colnum);
    out.lon_tangent = nan(1,colnum);
    out.lat_tangent  = nan(1,colnum);
    out.pressure_hPa = nan(rownum,colnum);
%     out.lon = nan(rownum,colnum);
%     out.lat = nan(rownum,colnum);
    out.quality_flags = nan(rownum,colnum);
    itotal = 0;
    for n = 1:length(datafolders) % go through the folders corresponding to each year and month       
        %% get a list of the .asc files in the directory
        tempdir = dir(strcat(datadirectory,'/',datafolders{n},'/','*uo3g*.dat')); % the names of the data files
        fileall = {tempdir.name};
        fprintf('%s\n',datafolders{n});
        fprintf('Reading files...\n')
        for i = 1:length(fileall)
            itotal = itotal + 1;
            filei = strcat(datadirectory,'/',datafolders{n},'/',fileall{i});
%             fprintf('%s\n',filei);
% % %             try % some files are blank
                [ maestroi ] = read_maestro_ascdata_v1p2(filei);
                if isstruct(maestroi)
                out.source_file{itotal} = maestroi.source_file;
                out.version{itotal} = maestroi.version;
                out.occultation(itotal) = maestroi.occultation;
                out.sr1ss0(itotal) = maestroi.sr1ss0;
                out.date_mjd(itotal) = maestroi.date_mjd;
                out.vmr(:,itotal) = maestroi.vmr;
                out.vmr_error(:,itotal) = maestroi.vmr_error;
                out.altitude_km(:,itotal) = maestroi.altitude_km; % all the altitude grids are the same
                out.quality_flags(:,itotal) = maestroi.quality_flags;
                else
                    warning('something went wrong with reading %s', filei)
                end
% % %             catch
% % %                 warning('something went wrong with reading %s', filei)
% % %             end
        end
    end
    out.source_file = out.source_file(1:itotal);
    out.version = out.version(1:itotal);
    tanstruct = reduce_tanstruct_by_rowindex(out, 1:itotal); % remove the columns that werent filled with data
    
    savedest = strcat(matdirectory,'/',savename);
    fprintf('\nsaving DMP data to %s\n', savedest);
    save(savedest,'tanstruct');
    fprintf('done\n')
else
    fprintf('\nThe directory, %s, does not exist. Please create it first\n', matdirectory)
end
end

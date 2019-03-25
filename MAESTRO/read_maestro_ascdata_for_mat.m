function [ tanstruct ] = read_maestro_ascdata_for_mat( path_to_maestro_data )
%A function to read the ACE MAESTRO O3 files and make a .mat
%structure. This function assumes that all of the MAESTRO files are in the
%same directory. The output folder for the .mat data should be specified by
%the user below. 

% *INPUT*
%           path_to_maestro_data: STRING - the path to the directory that
%           holds the data files. 
%
% *OUTPUT*
%           .mat files will be written to the outut folder that is
%           specified below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2018

%% Define some things
datadirectory = path_to_maestro_data;
if ~isdir(path_to_maestro_data)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',path_to_maestro_data)
    error('The directory containing the maestro files couldn''t be found')
end

%%USER DEFINED
% matdirectory = '/Users/niall/Dropbox/climatology/nryan/dmptest10'; % edit this to your output directory
matdirectory = 'C:\Users\ryann\ACE\MAESTRO\matdata\'; % edit this to your output directory
% matdirectory = 'F:\ACE\dmp3.5';

%%STANDARD
savename = 'MAESTRO_v3p13_O3';

%%
if isdir(matdirectory) 
    
    %% get a list of the folders in the directory
    tempdir = dir(datadirectory);
    datafolders = {tempdir.name}; % asuming that all that is in the directory are the folders with the data
    datafolders = datafolders(3:end); % get rid of ',' and'..' 
    
    %% Read out the data from each of the files and put them into a structure
    % There seems to be 311 levels of O3 data for now.
    rownum = 311;
    colnum = 2e5; % i think there's about 1.2e5 measurements from maestro at the moment, so this is an overshoot
    out.source_file = cell(1,colnum); %string
    out.version = cell(1,colnum); % string
%     out.source_file = cell(1,1);
%     out.version = cell(1,1);
    out.occultation = nan(1,colnum);
    out.sr1ss0 = nan(1,colnum);
    out.beta_angle = nan(1,colnum);
    out.date_mjd = nan(1,colnum);
    out.gas ='O3'; % this is fixed for now as I'm only doing O3 at the moment
    out.altitude_km = nan(rownum,colnum);
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
        tempdir = dir(strcat(datadirectory,'/',datafolders{n},'/','*vo3g*.dat')); % the names of the data files
        fileall = {tempdir.name};
        fprintf('%s\n',datafolders{n});
        fprintf('Reading files...\n')
        for i = 1:length(fileall)
            itotal = itotal + 1;
            filei = strcat(datadirectory,'/',datafolders{n},'/',fileall{i});
            %             fprintf('%s\n',filei);
            [ maestroi ] = read_maestro_ascdata(filei);
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
        end
    end
    out.source_file = out.source_file(1:itotal);
    out.version = out.version(1:itotal);
    tanstruct = reduce_tanstruct_by_rowindex(out, 1:itotal); % remove the columns that werent filled with data
    
    savedest = strcat(matdirectory,'/',savename);
    fprintf('\nsaving MAESTRO data to %s\n', savedest);
    save(savedest,'tanstruct');
    fprintf('done\n')
else
    fprintf('\nThe directory, %s, does not exist. Please create it first\n', matdirectory)
end
end

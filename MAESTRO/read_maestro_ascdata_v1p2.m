function [ maestro_out ] = read_maestro_ascdata_v1p2( path_to_file)
%A function to read the ACE MAESTRO v3.13 .asc, DMP data, and output the relevant
%information in a matlab structure. 

% *INPUT*    
%           filename: STRING - the name of the .txt file that holds the ACE
%           MAESTRO data
%
% *OUTPUT*
%           maestro_out: MATLAB STRUCTURE with th MAESTRO data. The fields
%           of the structure are the same as those as are used for the ACE
%           data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 05/2018 

%% Define some things
filein = path_to_file;
[~, fname, fext] = fileparts(filein); % for cases when the path to a file is also entered, as opposed to just the name of the file
filein = strcat(fname,fext);
out.source_file = filein;
out.gas = 'O3'; % this is fixed for now. I'm only doing O3 at the moment

% if nargin < 2
%     orbit_table_path = 'C:\Users\ryann\ACE\MAESTRO\SunriseSunsetTable';
%     temp = read_maestro_orbit_table(orbit_table_path);
%     orbit_info = temp.sr1ss0_table;
%     clear temp
% else
%     orbit_info = orbit_table;
% end

% dmpname = dir(filein); % gets the details of the .asc file in the directory
% dmpname = cell2mat({dmpname.name}); % gets the name of the .asc file in a string format
%% Read the filename info

% the occulation number and sunrise/sunset info from the name of the file
name_info = textscan(filein, '%s%s%s%s%s', 'Delimiter', '_'); % The file naming for MAESTRO data is 'ss2830_vo3g_040221_013059_27.dat'
orbit = cell2mat(name_info{1});
out.occultation = single(str2double(orbit(3:end))); % the occulation name is of the type: sr1234 or ss12345, for example
sr1ss0 = orbit(1:2); % get sunrise (sr) or sunset (ss) info. sr := 1, ss := 0.
if strcmp(sr1ss0,'sr')
    if strcmp(orbit(end),'a')
        out.sr1ss0 = 3;
    else
        out.sr1ss0 = 1;
    end
elseif strcmp(sr1ss0, 'ss')
    if strcmp(orbit(end),'a')
        out.sr1ss0 = 2;
    else
        out.sr1ss0 = 0;
    end
else
    disp('something is wrong with reading the occulation number')
end
out.beta_angle = nan; % don't have this info yet

% get the date and change into MJD
yearanddate = cell2mat(name_info{3});
year = 2000 + str2double(yearanddate(1:2));
month = str2double(yearanddate(3:4));
day = str2double(yearanddate(5:6));

time = cell2mat(name_info{4});
hour = str2double(time(1:2));
min = str2double(time(3:4));
sec = str2double(time(5:6));

out.date_mjd = date2mjd(year,month,day,hour,min,sec); % calculate the MJD

%% get the data from inside the file
% the data is in columns of 1:index  2:height  3:VMR  4:VMR_ERROR
% 5:retreived
fid1 = fopen(path_to_file);
for n = 1:6
    dummy = fgetl(fid1); %#ok<NASGU> % there are 10 header lines before the data
end
for n = 7
    tline = fgetl(fid1);
    if tline == -1
        warning('this appears to be an empty or corrupted file')
        fclose(fid1);
        maestro_out = 0;
        return
    end
    C = textscan(tline, '%s');
    out.version = cell2mat(C{1,1}(4,1));
end
for n = 8:10
    dummy = fgetl(fid1); %#ok<NASGU> % there are 10 header lines before the data
end
i = 0;
for n = 11:211 %the number of lines that contain the data table.
    i = i+1;
    tline = fgetl(fid1);
    C = textscan(tline, '%s');
    C = C{1,1};
%     tdata = cell2mat(C);
    out.altitude_km(i,1) = str2double(C{2});
    if strcmp(C{3}, '-1.#JE+000') % this is for bad data
        out.vmr(i,1) = nan;
        out.vmr_error(i,1) = nan;
    else
        out.vmr(i,1) = str2double(C{3});
        out.vmr_error(i,1) = str2double(C{4}).*out.vmr(i,1); % the error is fractional. we want vmr.
    end 
    retrieved = str2double(C{5});
    
    if isnan(out.vmr(i,1))
        out.quality_flags(i,1) = 7; % the ACE-FTS flag for an error in the data
    else
        if retrieved == 0
            out.quality_flags(i,1) = 8; % the ACE-FTS flag for a priori data
        elseif retrieved == 1
            out.quality_flags(i,1) = 0; % the ACE-FTS flag for data with no issues
        end
    end
end
fclose(fid1);

%% sort the data by altitude
[out.altitude_km, Ialt] = sort(out.altitude_km);
out.vmr = out.vmr(Ialt);
out.vmr_error = out.vmr_error(Ialt);
out.quality_flags = out.quality_flags(Ialt);
out.lat_tangent = nan;
out.lon_tangent = nan;

maestro_out = out;
%
end


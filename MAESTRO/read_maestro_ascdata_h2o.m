function [ maestro_out ] = read_maestro_ascdata_h2o( path_to_file)
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
%   NJR - 11/2018 

%% Define some things
filein = path_to_file;
[~, fname, fext] = fileparts(filein); % for cases when the path to a file is also entered, as opposed to just the name of the file
filein = strcat(fname,fext);
out.source_file = filein;
out.gas = 'H2O'; % this is fixed for now.
lalt = 50; % choose a value that will be larger than the number of levels

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
name_info = textscan(filein, '%s%s%s%s%s%s', 'Delimiter', '_'); % The file naming for MAESTRO data is 'ss2830_vo3g_040221_013059_27.dat'
% test = name_info
% return
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
out.version = [];

%% get the data from inside the file
% the data is in columns of 1:height  2:VMR  3:VMR_ERROR
out.altitude_km = nan(lalt,1);
out.vmr = nan(lalt,1);
out.vmr_error = nan(lalt,1);
out.lat = nan(lalt,1);
out.lon = nan(lalt,1);
out.quality_flags = nan(lalt,1);
fid1 = fopen(path_to_file);
for n = 1
    dummy = fgetl(fid1); %#ok<NASGU> % header lines before the data
end
for n = 2
    tline = fgetl(fid1);
    if tline == -1
        warning('this appears to be an empty or corrupted file')
        fclose(fid1);
        maestro_out = 0;
        return
    end
    C = textscan(tline, '%s');
%     test = C
%     return
    out.lat(:) = str2double(cell2mat(C{1,1}(4,1)));
    out.lon(:) = str2double(cell2mat(C{1,1}(5,1)));
end
% the h2o data has varying altitude levels and different numbers of levels
i = 0;
while 1
    i = i+1;
    tline = fgetl(fid1);
    C = textscan(tline, '%s');
    C = C{1,1};
    if strcmp(C{1},char(26)), break; end % end the loop if you get the square symbol
    out.altitude_km(i,1) = str2double(C{1});
    out.vmr(i,1) = str2double(C{2});
    out.vmr_error(i,1) = str2double(C{3});
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


function [ dmp_out ] = read_ace_dmpv2( filename )
%A function to read the ACE v2.0, .asc, DMP data, and output the relevant
%information in a matlab structure. 

% *INPUT*    
%           filename: STRING - the name of the .txt file that holds the ace
%           DMP data
%
% *OUTPUT*
%           dmp_out: MATLAB STRUCTURE with the following fields of ACE DMP
%           data: sourcefile, occultation, sr1ss0, date_mjd, lon_tangent,
%           lat_tangent, version, altitude_km, T, pressure_hPa, lon, lat,
%           Theta, spv, and eql.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 10/2017 

%% Define some things
filein = filename;
out.source_file = filein;
% dmpname = dir(filein); % gets the details of the .asc file in the directory
% dmpname = cell2mat({dmpname.name}); % gets the name of the .asc file in a string format
%% Read the file info
fid1 = fopen(filein);

% the occulation number and sunrise/sunset info
for n = 1
    tline = fgetl(fid1);
end
orbit = textscan(tline, '%s');
orbit = cell2mat(orbit{1,1}(2,1));
out.occultation = single(str2double(orbit(7:end))); % the occulation name is of the type: ace.sr1234 or ace.ss12345, for example
sr1ss0 = orbit(5:6); % get sunrise (sr) or sunset (ss) info. sr := 1, ss := 0.
if strcmp(sr1ss0,'sr')
    out.sr1ss0 = 1;
elseif strcmp(sr1ss0, 'ss')
    out.sr1ss0 = 0;
else
    disp('something is wrong with reading the occulation number')
end

% get the date and change into MJD
for n = 2
    tline = fgetl(fid1); % date is on the second line
end
dateandtime = textscan(tline, '%s');
date = cell2mat(dateandtime{1,1}(2,1));
year = str2double(date(1:4)); % date is of the type yyyy-mm-dd
month = str2double(date(6:7));
day = str2double(date(9:10));
time = cell2mat(dateandtime{1,1}(3,1));% time is of the type hh:mm:ss.ss+00
hour = str2double(time(1:2));
min = str2double(time(4:5));
%sec = str2double(time(7:8));
if length(time) == 14
    sec = str2double(time(7:11));
elseif length(time) == 13 % sometimes the time is of the type hh:mm:ss.s+00
    sec = str2double(time(7:10));
elseif length(time) == 11 %sometimes the time is of the type hh:mm:ss+00
    sec = str2double(time(7:8));
end
out.date_mjd = date2mjd(year,month,day,hour,min,sec); % calculate the MJD
% get the lat and lon
for n = 3
    tline = fgetl(fid1); % lon is on the third line
end
lon = textscan(tline, '%s');
out.lon_tangent = str2double(cell2mat(lon{1,1}(2,1)));
for n = 4
    tline = fgetl(fid1); % lat is on the fourth line
end
lat = textscan(tline, '%s');
out.lat_tangent = str2double(cell2mat(lat{1,1}(2,1)));

%discard the next line. it is the source information for the DMP: ACE GLC
%data
for n = 5
    dummy = fgetl(fid1); %#ok<NASGU> % dummy fith line
end

%get the model info that was used or the DMP
for n = 6
    tline = fgetl(fid1); % sixth line
end
version = textscan(tline, '%s');
out.version = strcat(cell2mat(version{1,1}(2,1)),cell2mat(version{1,1}(3,1)));

for n = 7:9
    dummy = fgetl(fid1); %#ok<NASGU> % lines seven to eleven aren't needed
end

for n = 10
    tline = fgetl(fid1); % tenth line
end
tropopauses = textscan(tline, '%s');
out.tropopauses_km = [str2double(tropopauses{1,1}(1,1)), str2double(tropopauses{1,1}(2,1)), str2double(tropopauses{1,1}(3,1)), str2double(tropopauses{1,1}(4,1))]';

for n = 11:12
    dummy = fgetl(fid1); %#ok<NASGU> % lines seven to eleven aren't needed
end


%% get the data
% %the columns contain information in the following order;
% % 1:z  2:ACE-T  3:ACE-p  4:ACE-lon  5:ACE-lat  6:Met-T  7:Met-Theta
% % 8:delh-T  9:GPH  10:Zon-Wind  11:Mer-Wind  12:PV  13:sPV  14:EqL  
% % 15:delh-PV  16:RelVor  17:DTDZ
% % There seems to be about 10 data points per day, on average, up to now. so
% % 20 years of data would be 73000 data points. make this an overestimate so
% % to have a size for the predefined matrices below, with 77 levels of dmp
% % data for now
% out.altitude_km = nan(77,73000); % see above
% out.T = nan(77,73000);
% out.pressure_hPa = nan(77,73000);
% out.lon = nan(76,73000);
% out.lat = nan(76,73000);
% out.theta = nan(76,73000);
% out.spv = nan(76,73000);
% out.eql = nan(76,73000);     ***********for the read_dmp_for_mat

% the columns contain information in the following order;
% 1:z  2:ACE-T  3:ACE-p  4:ACE-lon  5:ACE-lat  6:Met-T  7:Met-Theta
% 8:delh-T  9:GPH  10:Zon-Wind  11:Mer-Wind  12:PV  13:sPV  14:EqL  
% 15:delh-PV  16:RelVor  17:DTDZ
ndmp = 76; % 76 levels of dmp data for now
out.altitude_km = nan(ndmp,1); % see above
out.T = nan(ndmp,1);
out.pressure_hPa = nan(ndmp,1);
out.lon = nan(ndmp,1);
out.lat = nan(ndmp,1);
out.Theta = nan(ndmp,1);
out.spv = nan(ndmp,1);
out.eql = nan(ndmp,1);
i = 1;
for n = 13:88 %the length of the data columns: 76
    tline = fgetl(fid1);
    C = textscan(tline, '%f');
    dmpi = cell2mat(C);
    %dmpi = dmpi';
    out.altitude_km(i) = dmpi(1);
    out.T(i) = dmpi(2);
    out.pressure_hPa(i) = dmpi(3);
    out.lon(i) = dmpi(4);
    out.lat(i) = dmpi(5);
    out.Theta(i) = dmpi(7);
    out.spv(i) = dmpi(13);
    out.eql(i) = dmpi(14);
    i = i + 1;
end

fclose(fid1);
dmp_out = out;

end


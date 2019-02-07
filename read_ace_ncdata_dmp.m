function [ dmp_out ] = read_ace_ncdata_dmp( filename_glc, filename_trop )
%A function to read the ACE netcdf DMP data, and output the relevant
%information in a matlab structure. These DMPS were first prodiced in 2018
%by Luis Millan and Gloria Manney.

% *INPUT*    
%           filename_glc: STRING - the name of the .nc4 file that holds the
%           ace DMP data
%
%           filename_trop: STRING - the name of the .nc4 file that holds
%           the ace tropopause data.
%
% *OUTPUT*
%           dmp_out: MATLAB STRUCTURE with the following fields of ACE DMP
%           data: sourcefile, occultation, sr1ss0, date_mjd, lon_tangent,
%           lat_tangent, version, altitude_km, T, pressure_hPa, lon, lat,
%           Theta, spv, and eql.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 10/2018 

%% Define some things
filein = filename_glc;
out.source_file = filein;
% dmpname = dir(filein); % gets the details of the .asc file in the directory
% dmpname = cell2mat({dmpname.name}); % gets the name of the .asc file in a string format
%% Read the file info
% get the date and change into MJD
yearanddate = h5read(filein,'/Date'); % time x 1 cell
norbit = length(yearanddate);
hourstring = h5read(filein,'/Hour'); % time x 1 cell
out.date_mjd = nan(1,norbit);
rounddate_glc = nan(1,norbit);
for i = 1:norbit
    year = str2double(yearanddate{i}(1:4));
    month = str2double(yearanddate{i}(5:6));
    day = str2double(yearanddate{i}(7:8));
    hour = str2double(hourstring{i});
    out.date_mjd(i) = date2mjd(year, month, day, hour); % calculate the MJD
    rounddate_glc(i) = date2mjd(year, month, day, round2(hour,1e-4)); % calculate the MJD
end
out.version = cell2mat(h5read(filein, '/Met_Info'));
out.altitude_km = double(h5read(filein, '/Altitude'));
out.T = double(h5read(filein, '/Temperature'));
out.pressure_hPa = double(h5read(filein, '/Pressure'));
out.lon = double(h5read(filein, '/Lon'));
out.lat = double(h5read(filein, '/Lat'));
out.theta = double(h5read(filein, '/Theta'));
out.spv = double(h5read(filein, '/sPV'));
out.eql = double(h5read(filein, '/EqL'));

%% if the tropopause file is also input
if nargin > 1
    filein = filename_trop;
    %% make sure that the dates match the glc file
    % get the date and change into MJD
    yearanddate = h5read(filein,'/Date'); % time x 1 cell
    norbit = length(yearanddate);
    hour = double(h5read(filein,'/Hour')); % time x 1 cell
    trop.date_mjd = nan(1,norbit);
    rounddate_trop = nan(1,norbit);
    for i = 1:norbit
        year = str2double(yearanddate{i}(1:4));
        month = str2double(yearanddate{i}(5:6));
        day = str2double(yearanddate{i}(7:8));
        trop.date_mjd(i) = date2mjd(year, month, day, hour(i)); % calculate the MJD
        rounddate_trop(i) = date2mjd(year, month, day, round2(hour(i),1e-4)); % calculate the MJD
    end
    rounddate_glc = round2(rounddate_glc, 1e-4); % you need this because the glc file has the hours in text with 6 sig' digits
    rounddate_trop = round2(rounddate_trop, 1e-4);
%     vdate_glc = datevec(mjd2datenum(out.date_mjd));
%     vdate_trop = datevec(mjd2datenum(trop.date_mjd));
    if isequal(rounddate_glc, rounddate_trop) % if the dates match
        out.source_file_trop = filein;
        out.tropopauses_hPa_WMO = double(h5read(filein, '/Press_at_WMO_Tropopauses'));
        out.tropopauses_km_WMO = double(h5read(filein, '/Altitude_at_WMO_Tropopauses'));
        out.tropopauses_hPa_Dyn = double(h5read(filein, '/Press_at_Dyn_Tropopauses'));
        out.tropopauses_km_Dyn = double(h5read(filein, '/Altitude_at_Dyn_Tropopauses'));
    else
        disp('the dates/times of the GLC file and the TROP file don''t match. Not including the TROP data here')
    end
else
    out.tropopauses_hPa_WMO = nan(4,norbit);
    out.tropopauses_km_WMO = nan(4,norbit);
    out.tropopauses_hPa_Dyn = nan(5,4,norbit);
    out.tropopauses_km_Dyn = nan(5,4,norbit);
end

dmp_out = out;
%
end


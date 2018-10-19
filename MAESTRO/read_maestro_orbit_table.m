function [ data_out ] = read_maestro_orbit_table( filename )
%A function to read the ACE GLC files and put the data in a structure.
%*************** NOTE ***************
%This function should be rendered obsolete as there is a netcdf GLC file
%available now. The function 'merge_ace_glc.m' should be used instead.
%-NJR 11/17

% *INPUT*    
%           filename: STRING - the name of the .txt file that holds the ace
%           DMP data
%
% *OUTPUT*
%           data_out: MATLAB STRUCTURE with the following fields of ACE DMP
%           data: sourcefile, occultation, sr1ss0, date_mjd, lon_tangent,
%           lat_tangent, version, altitude_km, T, pressure_hPa, lon, lat,
%           Theta, spv, and eql.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 11/2017 

%% Define some things
filein = filename;
out.source_file = filein;

%% Read the file info
fid1 = fopen(filein);

for n = 1:3
    fgetl(fid1); % the first 3 lines are header info
end

%get the orbit info. the data is in columns of 'orbit_number date time latitude longitude beta-angle'.
i = 0;
while 1
    tline = fgetl(fid1);
    if ~ischar(tline)
        disp('End of list')
        disp(' ');
        break
    end
    i = i + 1;
    data_i = textscan(tline,'%s');
    
    %orbit
    orbit = cell2mat(data_i{1,1}(1,1));
    out.occultation(1,i) = single(str2double(orbit(3:end))); % the occulation name is of the type: sr1234 or ss12345, for example
    sr1ss0 = orbit(1:2); % get sunrise (sr) or sunset (ss) info. sr := 1, ss := 0.
    if strcmp(sr1ss0,'sr')
        if strcmp(orbit(end),'a')
            out.sr1ss0(1,i) = 3;
        else
            out.sr1ss0(1,i) = 1;
        end
    elseif strcmp(sr1ss0, 'ss')
        if strcmp(orbit(end),'a')
            out.sr1ss0(1,i) = 2;
        else
            out.sr1ss0(1,i) = 0;
        end
    else
        disp('something is wrong with reading the occulation number')
    end
    
    %date and time
    yearanddate = cell2mat(data_i{1,1}(2,1)); % format is '2003-10-01'
    year = str2double(yearanddate(1:4));
    month = str2double(yearanddate(6:7));
    day = str2double(yearanddate(9:10));
    
    time = cell2mat(data_i{1,1}(3,1)); % format is '00:16:25+00'
    hour = str2double(time(1:2));
    min = str2double(time(4:5));
    sec = str2double(time(7:8));
    
    out.date_mjd(1,i) = date2mjd(year,month,day,hour,min,sec); % calculate the MJD
    
    %latitude and longitude
    out.lat_tangent(1,i) = str2double(data_i{1,1}(4,1));
    out.lon_tangent(1,i) = str2double(data_i{1,1}(5,1));
    
    % beta angle
    out.beta_angle(1,i) = str2double(data_i{1,1}(6,1));
    
end

fclose(fid1);
data_out = out;

end


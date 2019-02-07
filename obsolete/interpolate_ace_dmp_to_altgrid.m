function [ tanstruct_out ] = interpolate_ace_dmp_to_altgrid( dmpstruct_in, altitude_grid )
%A function to interpolate the variables of an ACE DMP strcture to a
%user-defined pressure grid. The interpolation is performed in log-pressure
%space.
%
% *INPUT*
%           dmpstruct_in: STRUCTURE - containins the ACE dmp data. It is
%           usually created using 'read_ace_dmp' or 'read_ace_dmp_for_mat'. 
%            
%           pressure_grid: VECTOR - a vector containing the pressure grid
%           onto which you want to interpolate the ace and dmp variables.
%           This pressure should be in units of hPa.
%
% *OUTPUT*
%           dmpstruct_out: STRUCTURE - output has almost the same
%           fields as the input, but with the data interpolated onto pgrid.
%           The new pressure grid is exchanged for the old one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% define some things
interptype = 'pchip';
req_points = 4; % the number of points required for an interpolation
dmp = remove999_ace_dmp(dmpstruct_in); % removes the 999 values from the ace dmps
%dmp = dmpstruct_in;
lorbit = length(dmp.occultation);% this is the same as the length of p
zdmp = dmp.altitude_km; % this is the same altitude grid for all ace dmps
zgrid = altitude_grid;
if isvector(zgrid)
   if ~iscolumn(zgrid)
      error('The interpolation grid should be an array or a column vector') 
   end
end
lgrid = length(zgrid(:,1));
colgrid = length(zgrid(1,:));
if colgrid > 1 && colgrid ~= lorbit
    error('the number of interpolation grids does not match the number of occultations')
end

% the fields of the new data file that are the same
dmpout.source_file = dmp.source_file;
dmpout.occultation = dmp.occultation;
dmpout.sr1ss0 = dmp.sr1ss0;
dmpout.date_mjd = dmp.date_mjd;
dmpout.lon_tangent = dmp.lon_tangent;
dmpout.lat_tangent = dmp.lat_tangent;
dmpout.version = dmp.version;
% the fields of the new tanstruct file that will be interpolated
dmpout.T = nan(lgrid,lorbit);
dmpout.lon = nan(lgrid,lorbit);
dmpout.lat = nan(lgrid,lorbit);
dmpout.Theta = nan(lgrid,lorbit);
dmpout.spv = nan(lgrid,lorbit);
dmpout.eql = nan(lgrid,lorbit);
dmpout.pressure_hPa = nan(lgrid,lorbit);
dmpout.altitude_km = zgrid;
if length(dmp.pressure_hPa(1,:)) == 1 || length(dmp.altitude_km(1,:)) > 1
    warning('It seems like this structure has already been interpolated');
end

%% Interpolate the fields of the ace structure
fprintf('\ninterpolating the dmp data\n')
zgridi = zgrid(:,1); % set this initially outside the loop so it doesn't have to reset each time if pgrid is a vector
zdmpi = zdmp(:,1);
pdmpi = dmp.pressure_hPa(:,1);
logpdmpi = log(pdmpi);
for i = 1:lorbit
    if length(zgrid(1,:)) > 1 % If there are multiple grids to which to interpolate, i.e. pgrid is an array
       zgridi = zgrid(:,i); % this is only changed if pgrid is not a vector
    end
    if length(zdmp(1,:)) > 1 % If the structure already has multiple altitude grids for some reason i.e. zdmp is an array
       zdmpi = zdmp(:,i); % this is only changed if zace is not a vector
    end
    lgood = ~isnan(dmp.T(:,i)); % the indices of the T field that do not contain nans. nans will be placed where there is no data using 'remove999_ace_dmp.m'
    if sum(lgood) >= req_points % need at least 4 points to perform the interpolation
        Ti = dmp.T(lgood,i); % reduce the size of the vector
        zdmpi_s = zdmpi(lgood); % reduce the size as well
        %interpolate the fields in log-pressure space
        dmpout.T(:,i) = interp1(zdmpi_s,Ti,zgridi,interptype,nan);
    end
    %Repeat the above process for each variable. They may contain different
    %amounts of NaNs
    lgood = ~isnan(dmp.lon(:,i));
    if sum(lgood) >= req_points
        loni = dmp.lon(lgood,i);
        zdmpi_s = zdmpi(lgood);
        dmpout.lon(:,i) = interp1(zdmpi_s,loni,zgridi,interptype,nan);
    end
    lgood = ~isnan(dmp.lat(:,i));
    if sum(lgood) >= req_points
        lati = dmp.lat(lgood,i);
        zdmpi_s = zdmpi(lgood);
        dmpout.lat(:,i) = interp1(zdmpi_s,lati,zgridi,interptype,nan);
    end
    lgood = ~isnan(dmp.Theta(:,i));
    if sum(lgood) >= req_points
        Thetai = dmp.Theta(lgood,i);
        zdmpi_s = zdmpi(lgood);
        dmpout.Theta(:,i) = interp1(zdmpi_s,Thetai,zgridi,interptype,nan);
    end
    lgood = ~isnan(dmp.spv(:,i));
    if sum(lgood) >= req_points
        spvi = dmp.spv(lgood,i);
        zdmpi_s = zdmpi(lgood);
        dmpout.spv(:,i) = interp1(zdmpi_s,spvi,zgridi,interptype,nan);
    end
    lgood = ~isnan(dmp.eql(:,i));
    if sum(lgood) >= req_points
        eqli = dmp.eql(lgood,i);
        zdmpi_s = zdmpi(lgood);
        dmpout.eql(:,i) = interp1(zdmpi_s,eqli,zgridi,interptype,nan);
    end
    % the pressure is a little bit different
    if length(dmp.pressure_hPa(1,:)) > 1 % this should really be an array if you are interpolating in altitude space. You shouldn't be interpolating already interpolated data
        logpdmpi = log(dmp.pressure_hPa(:,i)); % this is only changed if it is an array
    end
    lgood = ~isnan(logpdmpi); % the indices of the P field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
    if sum(lgood) >= req_points % need at least 4 points to perform the interpolation
        logpdmpi = logpdmpi(lgood); % reduce the size of the vector
        zdmpi_s = zdmpi(lgood); % reduce the size as well
        logpdmpouti = interp1(zdmpi_s,logpdmpi,zgridi,interptype,nan);
        dmpout.pressure_hPa(:,i) = exp(logpdmpouti);
    end
end
fprintf('Done\n')
tanstruct_out = dmpout;
end


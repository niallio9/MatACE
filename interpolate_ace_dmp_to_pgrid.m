function [ tanstruct_out ] = interpolate_ace_dmp_to_pgrid( dmpstruct_in, pressure_grid )
%A function to interpolate the variables of an ACE DMP strcture to a
%user-defined pressure grid. The interpolation is performed in log-pressure
%space.
%
% *INPUT*
%           dmpstruct_in: STRUCTURE - containins the ACE dmp data. It is
%           usually created using 'read_ace_dmp' or 'read_ace_dmp_for_mat'. 
%            
%           pressure_grid: ARRAY - a vector or matrix containing the
%           pressure grid onto which you want to interpolate the ace and
%           dmp variables. This pressure should be in units of hPa.
%
% *OUTPUT*
%           dmpstruct_out: STRUCTURE - output has almost the same
%           fields as the input, but with the data interpolated onto pgrid.
%           The new pressure grid is exchanged for the old one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% define some things
interptype = 'pchip';
dmp = remove999_ace_dmp(dmpstruct_in);
%dmp = dmpstruct_in;
lorbit = length(dmp.occultation);% this is the same as the length of p
zdmp = dmp.altitude_km; % this is the same altitude grid for all ace dmps
pdmp = dmp.pressure_hPa;
% next line is to add to the ace pressure grid so that there are not values that
% are equal when doing the interpolation.
padd = cumsum(ones(size(zdmp)))*eps*1e6; % 1x150. order of 1e-10
pdmp = pdmp + padd;
logpdmp = log(pdmp);
pgrid = pressure_grid;
logpgrid = log(pgrid);
lgrid = length(pgrid(:,1));
colgrid = length(pgrid(1,:));
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
%include the interpolated altitiude grids
dmpout.altitude_km = nan(lgrid,lorbit);
% the fields of the new tanstruct file that will be interpolated
dmpout.T = nan(lgrid,lorbit);
dmpout.lon = nan(lgrid,lorbit);
dmpout.lat = nan(lgrid,lorbit);
dmpout.Theta = nan(lgrid,lorbit);
dmpout.spv = nan(lgrid,lorbit);
dmpout.eql = nan(lgrid,lorbit);
dmpout.pressure_hPa = pgrid;
if length(dmp.pressure_hPa(1,:)) == 1 || length(dmp.altitude_km(1,:)) > 1
    warning('It seems like this structure has already been interpolated');
end

%% Interpolate the fields of the ace structure
fprintf('\ninterpolating the dmp data\n')
logpgridi = logpgrid(:,1); % set this initially outside the loop so it doesn't have to reset each time if pgrid is a vector
zdmpi = dmp.altitude_km(:,1);
for i = 1:lorbit
    %i
    if length(logpgrid(1,:)) > 1 % If there are multiple grids to which to interpolate, i.e. pgrid is an array
       logpgridi = logpgrid(:,i); % this is only changed if pgrid is not a vector
    end
    lgood = ~isnan(dmp.T(:,i)); % the indices of the T field that do not contain nans. nans will be placed where there is no data using 'remove999_ace_dmp.m'
    if sum(lgood) > 4 % need at least 4 points to perform the interpolation
        Ti = dmp.T(lgood,i); % reduce the size of the vector
        logpdmpi = logpdmp(lgood,i); % reduce the size as well
        % the following two lines sets as nan the values of pgrid that
        % don't correspond to valid data to avoid warnings in the
        % interpolation function.
        logpgridi_s = logpgridi;
        logpgridi_s( logpgridi > logpdmpi(1) | logpgridi < logpdmpi(end) ) = nan;
        %interpolate the fields in log-pressure space
        dmpout.T(:,i) = interp1(logpdmpi,Ti,logpgridi_s,interptype,nan);
    end
    %Repeat the above process for each variable. They may contain different
    %amounts of NaNs
    lgood = ~isnan(dmp.lon(:,i));
    if sum(lgood) > 4
        loni = dmp.lon(lgood,i);
        logpdmpi = logpdmp(lgood,i);
        logpgridi_s = logpgridi;
        logpgridi_s( logpgridi > logpdmpi(1) | logpgridi < logpdmpi(end) ) = nan;
        dmpout.lon(:,i) = interp1(logpdmpi,loni,logpgridi_s,interptype,nan);
    end
    lgood = ~isnan(dmp.lat(:,i));
    if sum(lgood) > 4
        lati = dmp.lat(lgood,i);
        logpdmpi = logpdmp(lgood,i);
        logpgridi_s = logpgridi;
        logpgridi_s( logpgridi > logpdmpi(1) | logpgridi < logpdmpi(end) ) = nan;
        dmpout.lat(:,i) = interp1(logpdmpi,lati,logpgridi_s,interptype,nan);
    end
    lgood = ~isnan(dmp.Theta(:,i));
    if sum(lgood) > 4
        Thetai = dmp.Theta(lgood,i);
        logpdmpi = logpdmp(lgood,i);
        logpgridi_s = logpgridi;
        logpgridi_s( logpgridi > logpdmpi(1) | logpgridi < logpdmpi(end) ) = nan;
        dmpout.Theta(:,i) = interp1(logpdmpi,Thetai,logpgridi_s,interptype,nan);
    end
    lgood = ~isnan(dmp.spv(:,i));
    if sum(lgood) > 4
        spvi = dmp.spv(lgood,i);
        logpdmpi = logpdmp(lgood,i);
        logpgridi_s = logpgridi;
        logpgridi_s( logpgridi > logpdmpi(1) | logpgridi < logpdmpi(end) ) = nan;
        dmpout.spv(:,i) = interp1(logpdmpi,spvi,logpgridi_s,interptype,nan);
    end
    lgood = ~isnan(dmp.eql(:,i));
    if sum(lgood) > 4
        eqli = dmp.eql(lgood,i);
        logpdmpi = logpdmp(lgood,i);
        logpgridi_s = logpgridi;
        logpgridi_s( logpgridi > logpdmpi(1) | logpgridi < logpdmpi(end) ) = nan;
        dmpout.eql(:,i) = interp1(logpdmpi,eqli,logpgridi_s,interptype,nan);
    end
    if length(dmp.altitude_km(1,:)) > 1 % this should really only be a vector though. You shouldn't be interpolating already interpolated data
        zdmpi = dmp.altitude_km(:,i); % this is only changed if it is an array
    end
    lgood = ~isnan(zdmpi); % the indices of the P field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
    if sum(lgood) > 4 % need at least 4 points to perform the interpolation
        zdmpi = zdmpi(lgood); % reduce the size of the vector
        logpdmpi = logpdmp(lgood,i); % reduce the size as well
        logpgridi_s = logpgridi;
        logpgridi_s( logpgridi > logpdmpi(1) | logpgridi < logpdmpi(end) ) = nan;
        dmpout.altitude_km(:,i) = interp1(logpdmpi,zdmpi,logpgridi_s,interptype,nan);
    end
end
fprintf('Done\n')
tanstruct_out = dmpout;
end


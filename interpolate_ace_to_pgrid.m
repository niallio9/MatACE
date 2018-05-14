function [ tanstruct_out ] = interpolate_ace_to_pgrid( tanstruct_in, pressure_grid )
%A function to interpolate the variables of ACE an data structure to a
%user-defined pressure grid. The interpolation is performed in log-pressure
%space. It is recommended to apply the quality flags before using this
%function, with 'apply_ace_flags.m'. The interpolation will be weird if the
%-999 values are kept, so the flags should be applied anyway. The quality
%flags cannot be interpolated here.
%
% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%            
%           pressure_grid: ARRAY - a vector or matrix containing the
%           pressure grid onto which you want to interpolate the ace and
%           dmp variables. This pressure should be in units of hPa.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has almost the same
%           fields as the input, but with the data interpolated onto the
%           input pressure grid(s). The new pressure grid is included. The
%           quality flags are removed as they cannot be interpolated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 10/2017
% NJR - 02/2018

%% define some things
gas = remove_ace_surface_pressure(tanstruct_in); % there are often multiple surface pressure values, which messes with interpolation
gas = filter_ace_pressure(gas); % remove profiles that have zero values for pressures in them
% gas = tanstruct_in;
interptype = 'pchip'; %the spline interpolation method is giving bad results near the boundaries
req_points = 4; % the number of points required for an interpolation
lorbit = length(gas.occultation);% this is the same as the length of p
zace = gas.altitude_km;
pace = gas.pressure_hPa;
%% Prepare the ACE pressure data so that it can used for interpolation
% next line is to add to the ace pressure grid so that there are not values that
% are equal when doing the interpolation.
padd = flipud(cumsum(ones(size(zace(:,1))))*eps*1e6); % 1x150. order of 1e-10
padd = repmat(padd,1,length(pace(1,:))); % needed to add this because matlab 2013b on deluge is causeing probelms when adding a vector to a matrix
pace = pace + padd;
clear padd
% padd5 = flipud(cumsum(ones(5,1))*eps*1e6); % 1x150. order of 1e-10. just adding to the lowest 5 levels to avoid multiple values of 1 
% pace(1:5,:) = pace(1:5,:) + padd5;
logpace = log(pace);
% padd = cumsum(ones(size(zace)))*eps*1e6; % 1x150. order of 1e-10
% logpace = logpace + padd;
pgrid = pressure_grid;
if isvector(pgrid)
   if ~iscolumn(pgrid)
      error('The interpolation grid should be an array or a column vector') 
   end
end
logpgrid = log(pgrid);
lgrid = length(pgrid(:,1));
colgrid = length(pgrid(1,:));
if colgrid > 1 && colgrid ~= lorbit
    error('The number of interpolation grids does not match the number of occultations')
end

% the fields of the new data file that are the same
gasout.source_file = gas.source_file;
gasout.occultation = gas.occultation;
gasout.sr1ss0 = gas.sr1ss0;
gasout.beta_angle = gas.beta_angle;
gasout.date_mjd = gas.date_mjd;
gasout.gas = gas.gas;
% the fields of the new tanstruct file that will be interpolated
gasout.vmr = nan(lgrid,lorbit);
% size(gasout.vmr)
gasout.vmr_error = nan(lgrid,lorbit);
%these two also stay the same
gasout.lon_tangent = gas.lon_tangent;
gasout.lat_tangent = gas.lat_tangent;
%include the interpolated altitiude grids
gasout.altitude_km = nan(lgrid,lorbit);
%add the new pressure grid
gasout.pressure_hPa = pgrid;
if isfield(gas,'pressure_hPa')
    if length(gas.pressure_hPa(1,:)) == 1 || length(gas.altitude_km(1,:)) > 1
        if length(gas.occultation) ~=1
            warning('It seems like this structure has already been interpolated');
        end
    end
end
if isfield(gas,'lat')
    gasout.lon = nan(lgrid,lorbit);
    gasout.lat = nan(lgrid,lorbit);
end
if isfield(gas,'eql')
   gasout.eql = nan(lgrid,lorbit); 
end
%% Interpolate the fields of the ace structure
fprintf('\nInterpolating the data...')
logpgridi = logpgrid(:,1); % set this initially outside the loop so it doesn't have to reset each time if pgrid is a vector
logpacei = logpace(:,1);
for i = 1:lorbit
    %i
    if mod(i,1000) == 0
        fprintf('\npast %i of %i\n', i, lorbit)
    end
    
    if length(logpgrid(1,:)) > 1 % If there are multiple grids to which to interpolate, i.e. pgrid is an array
       logpgridi = logpgrid(:,i); % this is only changed if pgrid is not a vector
    end
    if length(logpace(1,:)) > 1 % If the structure already has multiple altitude grids for some reason i.e. zace is an array
        logpacei = logpace(:,i); % this is only changed if zace is not a vector
    end
    lgood = ~isnan(gas.vmr(:,i)); % the indices of the vmr field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
    if sum(lgood) >= req_points % need at least 5 points to perform the interpolation
        vmri = gas.vmr(lgood,i); % reduce the size of the vector
        vmr_errori = gas.vmr_error(lgood,i);
        logpacei_s = logpacei(lgood); % reduce the size as well
        %interpolate the fields in log-pressure space
        gasout.vmr(:,i) = interp1(logpacei_s,vmri,logpgridi,interptype, nan);
        gasout.vmr_error(:,i) = interp1(logpacei_s,vmr_errori,logpgridi,interptype, nan);
    end
    if length(zace(1,:)) > 1 % this should really only be a vector though. You shouldn't be interpolating already interpolated data
        zacei = zace(:,i); %
    else
        zacei = zace;
        zacei(isnan(logpacei)) = nan; % this is now added because there are sometimes multiple surface pressure values (that are filtered out, above), and that information cannot be reflected in a single altitude profile
    end
    lgood = ~isnan(zacei); % the indices of the Z field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
    if sum(lgood) >= req_points % need at least 5 points to perform the interpolation
        zacei = zacei(lgood); % reduce the size of the vector
        logpacei_s = logpacei(lgood); % reduce the size as well
        %interpolate the fields in log-pressure space?????
        gasout.altitude_km(:,i) = interp1(logpacei_s,zacei,logpgridi, interptype, nan);
    end
    if isfield(gas,'lon')
        loni = gas.lon(:,i);
        loni(isnan(logpacei)) = nan;
        lgood = ~isnan(loni); % the indices of the vmr field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
        if sum(lgood) >= req_points % need at least 5 points to perform the interpolation
            loni = loni(lgood); % reduce the size of the vector
            logpacei_s = logpacei(lgood); % reduce the size as well
            %interpolate the fields in log-pressure space
            gasout.lon(:,i) = interp1(logpacei_s,loni,logpgridi,interptype, nan);
        end
        lati = gas.lat(:,i);
        lati(isnan(logpacei)) = nan;
        lgood = ~isnan(lati); % the indices of the vmr field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
        if sum(lgood) >= req_points % need at least 5 points to perform the interpolation
            lati = lati(lgood); % reduce the size of the vector
            logpacei_s = logpacei(lgood); % reduce the size as well
            %interpolate the fields in log-pressure space
            gasout.lat(:,i) = interp1(logpacei_s,lati,logpgridi,interptype, nan);
        end
    end
    if isfield(gas,'eql')
        eqli = gas.eql(:,i);
        eqli(isnan(logpacei)) = nan;
        lgood = ~isnan(eqli); % the indices of the vmr field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
        if sum(lgood) >= req_points % need at least 5 points to perform the interpolation
            eqli = eqli(lgood); % reduce the size of the vector
            logpacei_s = logpacei(lgood); % reduce the size as well
            %interpolate the fields in log-pressure space
            gasout.eql(:,i) = interp1(logpacei_s,eqli,logpgridi,interptype, nan);
        end
    end
end
disp('done')
tanstruct_out = gasout;

%
end


function [ tanstruct_latbin ] = bin_ace_by_lat( tanstruct_in, lat_bounds )
% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           lat_bounds: VECTOR - a vector containing the start and end
%           boundaries of the latitude bins. It is assumed that the end of
%           one bin corresponds to the start of another, measning that the
%           binning is performed for all data within the first and last
%           entry of the vector. For the standard ACE climatology,
%           lat_bounds will be -90:5:90.
%
% *OUTPUT*
%           tanstruct_latbin: STRUCTURE - output contains different
%           fields than the input. The fields relate to zonal mean values
%           of the input gas data. The output fields are listed below in
%           (A).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Check if the structure has latitude and longitude information
if ~isfield(tanstruct_in,'lat')
    error('there is no ''lat'' field in the input. You can add the lat/lon information with ''add_ace_glc.m'' ')
end

%% Define some things
gas = tanstruct_in;
latbnds = lat_bounds;
llat = length(latbnds) - 1;
if isfield(gas,'pressure_hPa') % from old version. no need to change 
    if isempty(gas.pressure_hPa)
        lalt = length(gas.altitude_km(:,1));
    else
        lalt = length(gas.pressure_hPa(:,1));
    end
end


%% Create a structure of with the required fields and size for the binned output
% this function also includes the grids from the ace data and chosen
% latitude bands performs some checks to make sure that the ace data and
% DMPs are for matching occultations, etc.. 
% the new fields are listed below in (A)
gaslatbin = make_ace_empty_latbinstruct( gas, latbnds );

%% Loop through the latitude bins and fill in the fields of the output structure
for i = 1:llat
    %i
    %get the reduced data for a given lat bin
    latstart = latbnds(i);
    latend = latbnds(i+1);
    warning off % supress warnings about reducing some structures to zero entries. will likely happen here sometimes
    [gas_lati] = subset_ace_by_lat(gas, latstart, latend); % this has the same fields as the input tanstruct.
    warning on
    % fill in the field arrays
    if ~isempty(gas_lati.occultation)
        gaslatbin.beta_angle_mean(i) = mean(gas_lati.beta_angle);
        gaslatbin.beta_angle_var(i) = std(gas_lati.beta_angle);
        gaslatbin.doy_mean(i) = mean(mjd2doy(gas_lati.date_mjd));
        gaslatbin.start_date = datestr(mjd2datenum(gas.date_mjd(1)));
        gaslatbin.end_date = datestr(mjd2datenum(gas.date_mjd(end)));
        gaslatbin.lon_tangent_mean(i) = nanmean(gas_lati.lon_tangent);
        gaslatbin.lat_tangent_mean(i) = nanmean(gas_lati.lat_tangent);
        gaslatbin.lon_mean(:,i) = nanmean(gas_lati.lon,2);
        gaslatbin.vmr_zonal(:,i) = nanmean(gas_lati.vmr,2);
        gaslatbin.vmr_zonal_var(:,i) = nanvar(gas_lati.vmr,[],2);
        gaslatbin.obs_count(:,i) = sum(~isnan(gas_lati.vmr), 2); % must calculate this here so you can use it for calculating errors next
        gaslatbin.vmr_zonal_error(:,i) = sqrt(nansum(gas_lati.vmr_error.^2, 2)) ./ gaslatbin.obs_count(:,i); % the propagated error on the VMRs
        gaslatbin.vmr_zonal_standard_error(:,i) = sqrt(gaslatbin.vmr_zonal_var(:,i) ./ gaslatbin.obs_count(:,i)); % the standard error ont he mean
        [gaslatbin.lst_median(i), gaslatbin.lst_mean(i), gaslatbin.lst_max(i), gaslatbin.lst_min(i), ...
            gaslatbin.lst_var(i)] = get_LST_stats(mjd2lst(gas_lati.date_mjd, gas_lati.lon_tangent));
        for j = 1:lalt
            %get the time lat and lon of each point that you use. the non
            %nan ones
%             j
%             ~isnan(gas_lati.vmr(j,:))
            obs_lon =  gas_lati.lon(j, ~isnan(gas_lati.vmr(j,:))); 
            obs_lat =  gas_lati.lat(j, ~isnan(gas_lati.vmr(j,:)));
            obs_date =  gas_lati.date_mjd(~isnan(gas_lati.vmr(j,:)));
            gaslatbin.obs_location{j,i} = [obs_lon;obs_lat;obs_date];
        end
    end
end

tanstruct_latbin = gaslatbin;
%
end

%% (A)
% fieldnames of the structure produced by 'bin_ace_by_lat'
% 
%     'source_file'
%     'beta_angle_mean'
%     'beta_angle_var'
%     'date_mjd_mean'
%     'gas'
%     'lon_tangent_mean'
%     'lat_tangent_mean'
%     'pressure_hPa'
%     'altitude_km_mean'
%     'vmr_zonal'
%     'vmr_zonal_var'
%     'vmr_zonal_error'
%     'vmr_zonal_standard_error'
%     'lat'
%     'lat_bounds'
%     'obs_count'
%     'obs_location'
%     'lst_median'
%     'lst_mean'
%     'lst_max'
%     'lst_min'
%     'lst_var'

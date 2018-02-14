function [ tanstruct_eqlbin ] = bin_ace_by_eql( tanstruct_in, eql_bounds )
% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           eql_bounds: VECTOR - a vector containing the start and end
%           boundaries of the equivalent latitude bins. It is assumed that
%           the end of one bin corresponds to the start of another, meaning
%           that the binning is performed for all data within the first and
%           last entry of the vector. For the standard ACE climatology,
%           eql_bounds will be -90:5:90.
%
% *OUTPUT*
%           tanstruct_eqlbin: STRUCTURE - output contains different
%           fields than the input. The fields relate to zonal mean values
%           of the input gas data. The output fields are listed below in
%           (A).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 02/2018

%% Check if the structure has latitude and longitude information
if ~isfield(tanstruct_in,'eql')
    error('there is no ''eql'' field in the input. You can add the eql information with ''merge_ace_dmp .m'' ')
end

%% Define some things
gas = tanstruct_in;
eqlbnds = eql_bounds;
llat = length(eqlbnds) - 1;
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
gaseqlbin = make_ace_empty_eqlbinstruct( gas, eqlbnds );

%% Loop through the latitude bins and fill in the fields of the output structure
for i = 1:llat
    %i
    %get the reduced data for a given lat bin
    eqlstart = eqlbnds(i);
    eqlend = eqlbnds(i+1);
    warning off % supress warnings about reducing some structures to zero entries. will likely happen here sometimes
    [gas_eqli] = subset_ace_by_eql(gas, eqlstart, eqlend); % this has the same fields as the input tanstruct.
    warning on
    % fill in the field arrays
    if ~isempty(gas_eqli.occultation)
        gaseqlbin.beta_angle_mean(i) = mean(gas_eqli.beta_angle);
        gaseqlbin.beta_angle_var(i) = std(gas_eqli.beta_angle);
        gaseqlbin.date_mjd_mean(i) = mean(gas_eqli.date_mjd);
        gaseqlbin.lon_tangent_mean(i) = mean(gas_eqli.lon_tangent);
        gaseqlbin.lat_tangent_mean(i) = mean(gas_eqli.lat_tangent);
        gaseqlbin.lon_mean(:,i) = nanmean(gas_eqli.lon,2);
        gaseqlbin.vmr_zonal(:,i) = nanmean(gas_eqli.vmr,2);
        gaseqlbin.vmr_zonal_var(:,i) = nanvar(gas_eqli.vmr,[],2);
        gaseqlbin.obs_count(:,i) = sum(~isnan(gas_eqli.vmr), 2); % must calculate this here so you can use it for calculating errors next
        gaseqlbin.vmr_zonal_error(:,i) = sqrt(nansum(gas_eqli.vmr_error.^2, 2)) ./ gaseqlbin.obs_count(:,i); % the propagated error on the VMRs
        gaseqlbin.vmr_zonal_standard_error(:,i) = sqrt(gaseqlbin.vmr_zonal_var(:,i) ./ gaseqlbin.obs_count(:,i)); % the standard error on the mean
        [gaseqlbin.lst_median(i), gaseqlbin.lst_mean(i), gaseqlbin.lst_max(i), gaseqlbin.lst_min(i), ...
            gaseqlbin.lst_var(i)] = get_LST_stats(mjd2lst(gas_eqli.date_mjd, gas_eqli.lon_tangent));
        for j = 1:lalt
            %get the time lat and lon of each point that you use. the non
            %nan ones
%             j
%             ~isnan(gas_lati.vmr(j,:))
            obs_lon =  gas_eqli.lon(j, ~isnan(gas_eqli.vmr(j,:))); 
            obs_eql =  gas_eqli.lat(j, ~isnan(gas_eqli.vmr(j,:)));
            obs_date =  gas_eqli.date_mjd(~isnan(gas_eqli.vmr(j,:)));
            gaseqlbin.obs_location{j,i} = [obs_lon;obs_eql;obs_date];
        end
    end
end

tanstruct_eqlbin = gaseqlbin;
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

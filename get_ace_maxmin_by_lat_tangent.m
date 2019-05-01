function [ maxvmr_latbins, minvmr_latbins, latbins ] = get_ace_maxmin_by_lat_tangent( tanstruct_in, lat_bounds, filter_in )
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
% NJR - 02/2019

%% Check if the structure has latitude and longitude information
% % if ~isfield(tanstruct_in,'lat')
% %     error('there is no ''lat'' field in the input. You can add the lat/lon information with ''add_ace_glc.m'' ')
% % end

%% Define some things
gas = tanstruct_in;
latbnds = lat_bounds;
llat = length(latbnds) - 1;
latbins = nan(1,llat);
for i = 1:llat
    latbins(i) = mean([latbnds(i),latbnds(i+1)]); %get the midpoints of the latitude bins
end
if isfield(gas,'pressure_hPa') % from old version. no need to change 
    if isempty(gas.pressure_hPa)
        lalt = length(gas.altitude_km(:,1));
    else
        lalt = length(gas.pressure_hPa(:,1));
    end
end
if nargin > 2
    filter = filter_in;
    if strcmp(filter,'fraction') == 1
        fraction_val = 0.0015;
        fraction_val = 0.01;
    end
else
    filter = 'none';
end


%% Create an array with the required fields and size for the output
maxvmr_latbins = nan(lalt, length(latbins));
minvmr_latbins = nan(lalt, length(latbins));

%% Loop through the latitude bins and fill in the max and min values by altitude
for i = 1:llat
    %i
    %get the reduced data for a given lat bin
    latstart = latbnds(i);
    latend = latbnds(i+1);
    warning off % supress warnings about reducing some structures to zero entries. will likely happen here sometimes
    [gas_lati] = subset_ace_by_lat_tangent(gas, latstart, latend); % this has the same fields as the input tanstruct.
    warning on
    % fill in the output arrays
    for j = 1: lalt
        vmr_lati_altj = gas_lati.vmr(j,:);
        vmr_lati_altj = vmr_lati_altj(~isnan(vmr_lati_altj)); % remove the nans
        if ~isempty(vmr_lati_altj)
            switch filter
                case 'fraction'
                    to_filter = round2(length(vmr_lati_altj) * fraction_val, 1); % get the number of values to filter out
                    vmr_lati_altj = sort(vmr_lati_altj);
                    vmr_lati_altj = vmr_lati_altj(to_filter : end - to_filter); % removes the values from the given vector
                otherwise
            end
            maxvmr_latbins(j, i) = max(vmr_lati_altj);
            minvmr_latbins(j, i) = min(vmr_lati_altj);
        end
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% might need this too: uses the the same (mostly largest) value for the whole of the resective poles
maxvmr_latbins(:, 1:5) = repmat(maxvmr_latbins(:,6), [1,5]);
maxvmr_latbins(:, end-4:end) = repmat(maxvmr_latbins(:,end-5), [1,5]);
minvmr_latbins(:, 1:5) = repmat(minvmr_latbins(:,6), [1,5]);
minvmr_latbins(:, end-4:end) = repmat(minvmr_latbins(:,end-5), [1,5]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% might need this one.
% maxvmr_latbins(:,1) = maxvmr_latbins(:,2);
% maxvmr_latbins(:,end) = maxvmr_latbins(:,end-1);
% minvmr_latbins(:,1) = minvmr_latbins(:,2);
% minvmr_latbins(:,end) = minvmr_latbins(:,end-1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
end
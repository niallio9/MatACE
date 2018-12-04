function [mlsstruct_acesample, chosen_rowcolumn] = sample_and_scale_mls_for_ace(mlsstruct_in, tanstruct_in, pratstruct, output_appendix)
%A function to sample MLS data according to the time/lat/lon/alt of ACE
%measurements. Data from a chemical box model is to scale the data.

% *INPUT*
%           mlsstruct_in: STRUCTURE - contains the gas specific MLS data.
%           Can be created with 'extract_mls_data.m'.
%
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           cmam_sample: STRUCTURE - with the data that has been sampled
%           from CMAM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 04/18
tic
%% Filter the ace data
ace = apply_ace_flags(tanstruct_in);
ace = filter_ace_pressure(ace);
ace = filter_ace_bad_lat(ace,10);% my standard choices for filtering the ACE data - NJR.
ace = filter_ace_bad_lon(ace,30);

%% Define some things
time_lim = 12/24; % the time colocation criteria, in fraction of day.
distance_lim = 1000; % the distance colocation criteria, in km.
Re = 6371; % radius of the earth
if ~isfield(tanstruct_in,'lat')
    error('There is no GLC lat/lon information in the ACE structure. Stopping');
end
mls = convert_mls_to_ace_format(mlsstruct_in); % do this so that you can use the ace functions on the mls data
% reduce ACE data to the date range of the MLS data
disp('subsetting ACE data to the MLS date range...')
Igood = find(ace.date_mjd >= min(mls.date_mjd) & ace.date_mjd <= max(mls.date_mjd));
ace = reduce_tanstruct_by_rowindex(ace, Igood);
size(ace.occultation)
disp('done')

nocc = length(ace.occultation);
lalt_ace = length(ace.altitude_km(:,1));
lalt_mls = length(mls.pressure_hPa(:,1));
prat = pratstruct; %#ok<NASGU>
pratmo_ratio_limit = Inf;
% pratmo_ratio_limit = pratmo_limit;
% vmr_limit = 2e-9;
% if nargin > 3
%     vmr_minmax = vmr_minmax_in;
% end
savedest = 'MLS_v4p2_ClO_acesample_12h_1000km_20042010';
if nargin > 3
    savedest = strcat(savedest,'_',output_appendix,'.mat');
else
    savedest = strcat(savedest,'.mat');
end
    

%% make a max and min for MLS data in latitude bounds
disp('finding min and max values for mls by latitude...')
latbnds = -90:5:90;
llatbins = length(latbnds) - 1;
latbins = nan(1,llatbins);
mls_minvmr_latbins = nan(lalt_mls,llatbins);
mls_maxvmr_latbins = nan(lalt_mls,llatbins);
for i = 1:llatbins
    latbins(i) = mean([latbnds(i),latbnds(i+1)]); %get the midpoints of the latitude bins
end
for i = 2 : llatbins - 1 % mls doesnt measure higher latitudes that 85 N/S
    mls_lati = subset_ace_by_lat_tangent(mls, latbnds(i), latbnds(i+1));
    mls_minvmr_latbins(:,i) = min(mls_lati.vmr,[],2); % by altitude
    mls_maxvmr_latbins(:,i) = max(mls_lati.vmr,[],2); % by altitude
end
clear mls_lati
mls_minvmr_latbins(:,1) = mls_minvmr_latbins(:,2); % mls doesnt measure higher latitudes that 85 N/S
mls_minvmr_latbins(:,end) = mls_minvmr_latbins(:,end-1);
mls_maxvmr_latbins(:,1) = mls_maxvmr_latbins(:,2); % mls doesnt measure higher latitudes that 85 N/S
mls_maxvmr_latbins(:,end) = mls_maxvmr_latbins(:,end-1);
disp('done')
%%
% mls_maxvmr_latbins = mls_maxvmr_latbins * 1.5; % adding this in to give some more freedom to the south pole
% mls_maxvmr_latbins = mls_maxvmr_latbins * inf;
%%
% test1 = mls_minvmr_latbins;
% test2 = mls_maxvmr_latbins;

%% output structure
out.source_file = mls.source_file;
out.occultation = ace.occultation;
out.sr1ss0 = ace.sr1ss0;
out.beta_angle = nan(1,nocc);
out.date_mjd = ace.date_mjd;
out.gas = mls.gas;
out.altitude_km = ace.altitude_km;
out.vmr = nan(size(ace.vmr));
out.vmr_error = nan(size(ace.vmr));
out.lat_tangent = ace.lat_tangent;
out.lon_tangent = ace.lon_tangent;
out.quality_flags = nan(size(ace.vmr));
out.pressure_hPa = ace.pressure_hPa;
out.lat = ace.lat;
out.lon = ace.lon;
out.distance = nan(size(ace.vmr));
out.time_diff = nan(size(ace.vmr));

chosen_rowcolumn = nan(lalt_ace.*nocc,2); % a matrix to hold the values of the chosen
% ratio_out = nan(lalt_ace, nocc);
out.lst_ratio = nan(size(ace.vmr));
k = 0;

%% go through ace occultations
disp('going through ACE occultations to find coincident MLS points...')
starton = 1;
for n = starton : nocc
    %         n
    if ~rem(n,10) || n == 1
        fprintf('past occultation %i of %i\n', n, nocc);
    end
    ace_n = reduce_tanstruct_by_rowindex(ace,n);
    lstace_n     = mjd2lst(ace_n.date_mjd, ace_n.lon); % vector of LSTs for nth ace occultation
    %% find the min and max vmr limits from the mls binned data
    [~,Imin] = min(ace_n.lat_tangent - latbins);
    if ace_n.lat_tangent >= latbnds(Imin)
        latbin_number = Imin;
    else
        latbin_number = Imin - 1;
    end
    vmr_limit_min = mls_minvmr_latbins(:,latbin_number);
    vmr_limit_max = mls_maxvmr_latbins(:,latbin_number);
    vmr_limit_min = interp1(mls.pressure_hPa(:,1), vmr_limit_min, ace.pressure_hPa(:,n)); %interpolate these limits to the ace pressure grid
    vmr_limit_max = interp1(mls.pressure_hPa(:,1), vmr_limit_max, ace.pressure_hPa(:,n));
    
    %% subset mls data to 'time_lim' hours around the ace data
    Itime = find(abs(mls.date_mjd - ace_n.date_mjd) < time_lim);
    if ~isempty(Itime)
        %                 fprintf('\n%i MLS profiles within %0.2f days\n', length(Itime), time_lim)
        mls3h = reduce_tanstruct_by_rowindex(mls, Itime); % % get the mls data that is within time_lim of the ace time. I3h size data structure
        %% Locate the biggest possible lat-lon range and then subset the data to that range. NJR - 07/2018 
        minlon_acen = min(ace_n.lon) * (pi/180);
        maxlon_acen = max(ace_n.lon) * (pi/180);
        minlat_acen = min(ace_n.lat) * (pi/180);
        maxlat_acen = max(ace_n.lat) * (pi/180);
        % use the distance criteria to find a lat/lon range
        dlat = distance_lim ./ Re; % from l = r * dlat. the one for longitude is: l = r*cos(lat) * dlon
        dlon = distance_lim ./ (Re * cos(min(abs([minlat_acen, maxlat_acen])))); % use the lat that is closet to the equator
        latrange = [minlat_acen - dlat, maxlat_acen + dlat] * (180/pi);
        lonrange = [minlon_acen - dlon, maxlon_acen + dlon] * (180/pi);
        mls3h = subset_ace_by_lat(mls3h, latrange(1), latrange(2)); % subset to latitude range
        mls3h = subset_ace_by_lon(mls3h, lonrange(1), lonrange(2)); % subset to the longitude range 
        
        if ~isempty(mls3h.occultation)
            %% interpolate all to the ace altitudes
            evalc('mls3h_int = interpolate_ace_to_pgrid(mls3h, ace.pressure_hPa(:,n));'); % interpolate the mls data to the pressure levels of the ace occultation. suppress command line output.
            mls3h_int.occultation(:) = ace_n.occultation;
            mls3h_int.sr1ss0(:) = ace_n.sr1ss0;
            mls3h_int.altitude_km = ace_n.altitude_km; % replace the interpolated altitude vector (contains nans) with the ace altitude vector. mls data was interpolated to this anyway, above.
            mls3h_int.quality_flags = nan(size(mls3h_int.vmr)); % add in a dummy quality flags matrix. The pratmo calculations later call 'apply_ace_flags.m'
            % get the distance between the ace measurement and the mls
            % measurements, at each altitude, in km.
            d = latlon2distance(mls3h_int.lat,mls3h_int.lon,ace_n.lat, ace_n.lon, mls3h_int.altitude_km); % ace altitude x I3h array
            Id = find(d <= distance_lim); % linear index here
            if ~isempty(Id)
                %                         fprintf('%i data points lie within %f km\n', length(Id), distance_lim)
                % subset to the data that is within 'distance_lim' of the ace measurement
                mls3h_int_d = reduce_tanstruct_data_by_index(mls3h_int, Id); % this will give an empty matrix if there are no points that fit the criteria.
                %% Get the chemical box model data for the occultation
                lstace_n_rep = repmat(lstace_n, size(mls3h_int_d.vmr(1,:))); %#ok<NASGU> % duplicate the lst vector to be the same number of profiles as mls measurements
                evalc('[~, ratio_prat3h_int_d] = scale_ace_with_pratmo_vmrs(mls3h_int_d, prat, lstace_n_rep);'); % get the pratmo ratios for the nth occultation and each of the mls profiles
                % we now need to recalculate the distances because the last reduction of the data above might
                %have removed some measurements.
                d = latlon2distance(mls3h_int_d.lat, mls3h_int_d.lon,ace_n.lat, ace_n.lon, mls3h_int_d.altitude_km); % ace altitude x I3h_d array
                % this new d will have all values less than the distance limit.
                % choose the closest values at each altitude. Have to go row by
                % row because the 'min' function can't deal with nans properly.
                for j = 1:length(mls3h_int_d.altitude_km(:,1))
                    % % % %             for j = 13:52 % these are the altitude levels at which there are no pressure values below 146 or above 1
                    %                                 j
                    % construct a new mls profile measurement from each of
                    % the closest measurements
                    dj = d(j,:); % the d values at the jth altitude. vector
                    ratio_prat3h_int_dj = ratio_prat3h_int_d(j,:); %#ok<NODEF> % pratmo values at the jth altitude
                    mls3h_int_dj = mls3h_int_d; % initialise a structure that will be altered into a subset of 'mls3h_int_d' in the loop below
                    %% this part is to deal with situations when MLS vmr is negative
                    ifound = find(mls3h_int_dj.vmr(j,:) < 0 & ratio_prat3h_int_dj > 1); % vmr < 0, pratmo > 1
                    if ~isempty(ifound)
%                         disp('making some negative vmrs positive for >1 scaling.')
%                         fprintf('n = %i. j = %i\n', n, j)
                        mls3h_int_dj.vmr(j,ifound) = mls3h_int_dj.vmr(j,ifound) * -1; % want to positively increase the vmr here, not make it a larger negative number
                    end
                    ifound = find(mls3h_int_dj.vmr(j,:) < 0 & ratio_prat3h_int_dj <= 1); % vmr < 0, pratmo <= 1
                    if ~isempty(ifound)
%                         disp('making some negative vmrs stay as they are.')
%                         fprintf('n = %i. j = %i\n', n, j)
                        mls3h_int_dj.vmr(j,ifound) = mls3h_int_dj.vmr(j,ifound) ./ ratio_prat3h_int_dj(ifound); % want to just keep the negative vmr the same here: assume there is no ClO and MLS has a negative measurement
                    end
                    
                    %% main loop for finding a value at (n,j)
                    vmr_x_prat_dj = mls3h_int_dj.vmr(j,:) .* ratio_prat3h_int_dj;
                    redo = 1; % to start the while loop
                    while nansum(dj) ~= 0 && nansum(ratio_prat3h_int_dj) ~= 0 && redo == 1
                        % reduce the values to those that fit the following
                        % crtiteria
                        inonan = find(~isnan(dj) & ~isnan(ratio_prat3h_int_dj) & ratio_prat3h_int_dj < pratmo_ratio_limit);% & vmr_x_prat_dj >= vmr_limit_min(j) & vmr_x_prat_dj <= vmr_limit_max(j)); % the indexes of the columns of dj that aren't nans
%                         size(dj)
                        dj = dj(inonan); % remove nan values
%                         size(dj)
                        if ~isempty(dj)
                            mls3h_int_dj = reduce_tanstruct_by_rowindex(mls3h_int_dj, inonan); % reduce the data to the rows that aren't nans at that altitude
                            ratio_prat3h_int_dj = ratio_prat3h_int_dj(inonan); % reduce the ratios too. 1 x length(dj)
                            vmr_x_prat_dj = vmr_x_prat_dj(inonan); % reduce the vmrs * ratios too. 1 x length(dj)
                            [~,idjmin] = min(dj,[],2); % get the index of the smallest distance for altitude j.
                            %check if the value has been chosen before
                            chosen_rowcolumn_n = [j, find(mls.date_mjd == mls3h_int_dj.date_mjd(idjmin))]; % the row column of the chosen point in the main data (interpolated)
                            already_chosen = ismember(chosen_rowcolumn_n, chosen_rowcolumn, 'rows'); % this will be [1,1] if the value has already been chosen
                            if already_chosen == 1
                                dj(idjmin) = nan; % change the value of dj to a nan and go back and choose another value
                                redo = 1;
                                disp('value already chosen, searching again')
                                fprintf('n = %i. j = %i\n', n, j)
                            else
                                redo = 0; % you have found a value that has not already been chosen so go ahead
                                % fill in the values for the chosen index
                                out.date_mjd_mls(j,n)     = mls3h_int_dj.date_mjd(idjmin);
                                out.vmr(j,n)        = mls3h_int_dj.vmr(j,idjmin);
                                out.vmr_error(j,n)  = mls3h_int_dj.vmr_error(j,idjmin);
                                %                                 out.lat_tangent(j,n)  = mls3h_int_dj.lat_tangent(idjmin);
                                %                                 out.lon_tangent(j,n)  = mls3h_int_dj.lon_tangent(idjmin);
                                out.lat_mls(j,n)        = mls3h_int_dj.lat(j,idjmin);
                                out.lon_mls(j,n)        = mls3h_int_dj.lon(j,idjmin);
                                out.distance(j,n)   = dj(idjmin);
                                out.time_diff(j,n) = out.date_mjd_mls(j,n) - out.date_mjd(n);
                                % record the rows and altitudes of the chosen value to
                                % make sure that it doesn't repeat
                                k = k+1;
                                chosen_rowcolumn(k,:) = chosen_rowcolumn_n;
                                out.lst_ratio(j,n) = ratio_prat3h_int_dj(idjmin);
                            end
                        end
                    end
                end
            else
                fprintf('no data points lie within %f km. n = %i.\n', distance_lim, n)
                % the nth 'out' column will remain all NaNs in this case
            end
        else
            fprintf('no data points lie within the calculated lat/lon range. n = %i.\n', n)
        end
    else
        fprintf('no data points lie within %f days. n = %i.\n', time_lim, n)
        % the nth 'out' column will remain all NaNs in this case
    end
    if mod(n,1000) == 0 || n == 1
%         partialsave = sprintf('mlsstruct_acesample_%i_%i', starton, n);
        mlsstruct_acesample_partialsave = out; %#ok<NASGU>
        disp('partial save made without scaling applied')
        save(savedest,'mlsstruct_acesample_partialsave');
    else
    end
end
out.vmr = out.vmr .* out.lst_ratio;
out.vmr_error = out.vmr_error .* out.lst_ratio;

mlsstruct_acesample = out;
fprintf('\nSaving sampled data with scaling applied to %s\n', savedest);
save(savedest,'mlsstruct_acesample')
disp('done')
%
toc
end

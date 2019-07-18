function [smilesstruct_acesample, chosen_rowcolumn] = sample_and_scale_mls_for_ace_no_vortex(mlsstruct_in, tanstruct_in, gas_name, pratstruct, output_appendix, minmax_mls)
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
%           pratstruct: STRUCTURE - contains the gas specific data for VMR
%           as a function of local solar time. Created using the script
%           'make_ace_gas_vmrs_with_pratmo(_by_year).mat'
%
%           output_appendix: STRING - an appendix to the name of the saved
%           output file. OPTIONAL.
%
%           minmax_mls: STRUCTURE - contains the gas specific information
%           about the maximum and minumum allowed values for gas at each
%           altitude and latitude bin. This is a rather specific structure.
%           The fields of the structure can be created using
%           'get_ace_maxmin_by_lat_tangent.m'. OPTIONAL
%
% *OUTPUT*
%           mlsstruct_acesample: STRUCTURE - with the data that has been
%           sampled from MLS.
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
% size(ace.occultation)
disp('done')

nocc = length(ace.occultation);
lalt_ace = length(ace.altitude_km(:,1));
lalt_mls = length(mls.pressure_hPa(:,1));

prat = pratstruct; %#ok<NASGU>
pratmo_ratio_limit = 10;
fprintf('A maximum scaling ratio of %f will be used\n', pratmo_ratio_limit)
% pratmo_ratio_limit = pratmo_limit;
% vmr_limit = 2e-9;
% if nargin > 3
%     vmr_minmax = vmr_minmax_in;
% end

savedest = sprintf('SMILES_v3p2_%s_acesample_12h_1000km', gas_name);
if nargin > 4
    if ~isempty(output_appendix)
        output_appendix = strcat('_',output_appendix);
    end
    savedest = strcat(savedest,output_appendix,'.mat');
else
    savedest = strcat(savedest,'.mat');
end


%% make a max and min for MLS data in latitude bounds
if nargin > 5
    if ~isempty(minmax_mls)
        if isscalar(minmax_mls)
            use_minmax = 1;
            lat_bounds = -90 : 5 : 90;
            lat_bins = -87.5 : 5 : 87.5;
            mls_maxvmr_latbins = ones(lalt_mls, length(lat_bins)) * minmax_mls;
%             size(mls_maxvmr_latbins)
            fprintf('hard maximum of %fppb chosen for the output\n', minmax_mls * 1e9)
        elseif isstruct(minmax_mls)
            use_minmax = 1;
            lat_bounds = minmax_mls.lat_bounds;
            lat_bins = minmax_mls.lat_bins;
            mls_maxvmr_latbins = minmax_mls.max_val;
            disp('hard maxima chosen for the output')
            disp('')
        else
            error('unknown style of input for minmax_mls')
        end
    else
        use_minmax = 0;
        disp('no hard maxima chosen for the output')
        disp('')
    end
else
    use_minmax = 0;
    disp('no hard maxima chosen for the output')
    disp('')
end

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
out.date_mjd_mls = nan(size(ace.vmr));
% out.time_diff = nan(size(ace.vmr));

chosen_rowcolumn = nan(lalt_ace.*nocc,2); % a matrix to hold the values of the chosen
% ratio_out = nan(lalt_ace, nocc);
out.lst_ratio = nan(size(ace.vmr));
k = 0;

%% go through ace occultations
disp('going through ACE occultations to find coincident MLS points...')
disp('')
starton = 1;
for n = starton : nocc
    %         n
    if ~rem(n,100) || n == 1
        fprintf('past occultation %i of %i\n', n, nocc);
    end
    ace_n = reduce_tanstruct_by_rowindex(ace,n);
    lstace_n     = mjd2lst(ace_n.date_mjd, ace_n.lon); % vector of LSTs for nth ace occultation
    
    if use_minmax == 1
        %% find the min and max vmr limits from the mls binned data
        [~,Imin] = min(ace_n.lat_tangent - lat_bins);
        if ace_n.lat_tangent >= lat_bounds(Imin)
            latbin_number = Imin;
        else
            latbin_number = Imin - 1;
        end
%         vmr_limit_min = mls_minvmr_latbins(:,latbin_number);
        vmr_limit_max = mls_maxvmr_latbins(:,latbin_number); % the variable is loaded above
%         vmr_limit_min = interp1(mls.pressure_hPa(:,1), vmr_limit_min, ace.pressure_hPa(:,n)); %interpolate these limits to the ace pressure grid
        vmr_limit_max = interp1(mls.pressure_hPa(:,1), vmr_limit_max, ace.pressure_hPa(:,n));
    else
        vmr_limit_max = Inf(size(ace.pressure_hPa(:,n)));
%         vmr_limit_min = -Inf;
    end
    
    %% subset mls data to 'time_lim' hours around the ace data
    Itime = find(abs(mls.date_mjd - ace_n.date_mjd) < time_lim);
    if ~isempty(Itime)
        %                 fprintf('\n%i MLS profiles within %0.2f days\n', length(Itime), time_lim)
        warning off
        mls_h = reduce_tanstruct_by_rowindex(mls, Itime); % % get the mls data that is within time_lim of the ace time. I_h size data structure
        warning on
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
        warning off
        mls_h = subset_ace_by_lat(mls_h, latrange(1), latrange(2)); % subset to latitude range
        mls_h = subset_ace_by_lon(mls_h, lonrange(1), lonrange(2)); % subset to the longitude range
        warning on
        
        if ~isempty(mls_h.occultation)
            %% interpolate all to the ace altitudes
            evalc('mls_h_int = interpolate_ace_to_pgrid(mls_h, ace.pressure_hPa(:,n));'); % interpolate the mls data to the pressure levels of the ace occultation. suppress command line output.
            clear mls_h % the last time is used for occultation n
            mls_h_int.occultation(:) = ace_n.occultation;
            mls_h_int.sr1ss0(:) = ace_n.sr1ss0;
            mls_h_int.altitude_km = ace_n.altitude_km; % replace the interpolated altitude vector (contains nans) with the ace altitude vector. mls data was interpolated to this anyway, above.
            mls_h_int.quality_flags = nan(size(mls_h_int.vmr)); % add in a dummy quality flags matrix. The pratmo calculations later call 'apply_ace_flags.m'
            %% get the distance between the ace measurement and the mls and subset to lat lon limits found above
            % measurements, at each altitude, in km.
            d = latlon2distance(mls_h_int.lat,mls_h_int.lon,ace_n.lat, ace_n.lon, mls_h_int.altitude_km); % ace altitude x I_h array
            Id = find(d <= distance_lim); % linear index here
            
            if ~isempty(Id)
% % %                 %                         fprintf('%i data points lie within %f km\n', length(Id), distance_lim)
% % %                 % subset to the data that is within 'distance_lim' of the ace measurement
                warning off
                mls_h_int_d_v = reduce_tanstruct_data_by_index(mls_h_int, Id); % this will give an empty matrix if there are no points that fit the criteria.
                warning on
% % %                 %% get the locations of the ace data with respect to a polar vortex and subset mls data to match
% % %                 vortex_loc_ace = get_ace_vortex_position(ace_n); % ace altitude x 1 vector
% % %                 vortex_loc_ace = repmat(vortex_loc_ace, size(mls_h_int_d.vmr(1,:))); % ace altitude x I_h_d array
% % %                 vortex_loc_mls = get_ace_vortex_position(mls_h_int_d);
% % %                 Iv = find(vortex_loc_ace - vortex_loc_mls == 0);
                
                if ~isempty(Id)
% % %                     warning off
% % %                     mls_h_int_d_v = reduce_tanstruct_data_by_index(mls_h_int_d, Id);
% % %                     warning on
                    %% Get the chemical box model data for the occultation
                    lstace_n_rep = repmat(lstace_n, size(mls_h_int_d_v.vmr(1,:))); %#ok<NASGU> % duplicate the lst vector to be the same number of profiles as mls measurements
                    evalc('[~, ratio_prat_h_int_d_v] = scale_ace_with_pratmo_vmrs(mls_h_int_d_v, prat, lstace_n_rep);'); % get the pratmo ratios for the nth occultation and each of the mls profiles
                    % we now need to recalculate the distances because the last reduction of the data above might
                    %have removed some measurements.
                    d = latlon2distance(mls_h_int_d_v.lat, mls_h_int_d_v.lon,ace_n.lat, ace_n.lon, mls_h_int_d_v.altitude_km); % ace altitude x I_h_d array
                    % this new d will have all values less than the distance limit.
                    % choose the closest values at each altitude. Have to go row by
                    % row because the 'min' function can't deal with nans properly.
                    
                    for j = 1:length(mls_h_int_d_v.altitude_km(:,1))
                        % construct a new mls profile measurement from each of
                        % the closest measurements
                        dj = d(j,:); % the d values at the jth altitude. vector
                        ratio_prat_h_int_d_vj = ratio_prat_h_int_d_v(j,:); %#ok<NODEF> % pratmo values at the jth altitude
                        mls_h_int_d_vj = mls_h_int_d_v; % initialise a structure that will be altered into a subset of 'mls_h_int_d_v' in the loop below
                        
                        %% this part is to deal with situations when MLS vmr is negative
                        % subbing in this new part to replace the two
                        % conditions below. Now what we want to do if vmr
                        % is negative and pratmo is >1, is to scale the vmr
                        % as if it were positive and subtract the initial
                        % vmr value. So scale the value into the positive,
                        % but keep the intial negative part too.
%                         size(mls_h_int_d_vj.vmr(j,:))
%                         whos
                        ifound = find(mls_h_int_d_vj.vmr(j,:) < 0 & ratio_prat_h_int_d_vj > 1); % vmr < 0, pratmo > 1
                        if ~isempty(ifound)
                            %                         disp('making some negative vmrs positive for >1 scaling.')
                            %                         fprintf('n = %i. j = %i\n', n, j)
                            mls_h_int_d_vj.vmr(j,ifound) = mls_h_int_d_vj.vmr(j,ifound) * -1; % change the vmr from negative to positive here. this is the same as x + 2*abs(x)
                            mls_h_int_d_vj.vmr(j,ifound) = ( (ratio_prat_h_int_d_vj(ifound) .* mls_h_int_d_vj.vmr(j,ifound)) - 2*mls_h_int_d_vj.vmr(j,ifound)) ./ ratio_prat_h_int_d_vj(ifound);
                            % this line above gives:
                            % ((pratmo*vmr) - 2*vmr) / pratmo
                            % which, when multiplied by pratmo later on,
                            % will give (pratmo*vmr) - 2*vmr)
                        end
                        
%                         ifound = find(mls_h_int_d_vj.vmr(j,:) < 0 & ratio_prat_h_int_d_vj > 1); % vmr < 0, pratmo > 1
%                         if ~isempty(ifound)
%                             %                         disp('making some negative vmrs positive for >1 scaling.')
%                             %                         fprintf('n = %i. j = %i\n', n, j)
%                             mls_h_int_d_vj.vmr(j,ifound) = mls_h_int_d_vj.vmr(j,ifound) * -1; % want to positively increase the vmr here, not make it a larger negative number
%                         end
%                         ifound = find(mls_h_int_d_vj.vmr(j,:) < 0 & ratio_prat_h_int_d_vj <= 1); % vmr < 0, pratmo <= 1
%                         if ~isempty(ifound)
%                             %                         disp('making some negative vmrs stay as they are.')
%                             %                         fprintf('n = %i. j = %i\n', n, j)
%                             mls_h_int_d_vj.vmr(j,ifound) = mls_h_int_d_vj.vmr(j,ifound) ./ ratio_prat_h_int_d_vj(ifound); % want to just keep the negative vmr the same here: assume there is no ClO and MLS has a negative measurement
%                         end
                        clear ifound
                        
                        %% main loop for finding a value at (n,j)
                        vmr_x_prat_dj = mls_h_int_d_vj.vmr(j,:) .* ratio_prat_h_int_d_vj;
                        redo = 1; % to start the while loop
                        while nansum(dj) ~= 0 && nansum(ratio_prat_h_int_d_vj) ~= 0 && redo == 1
                            % reduce the values to those that fit the following
                            % crtiteria
                            inonan = find(~isnan(dj) & ~isnan(ratio_prat_h_int_d_vj) & ratio_prat_h_int_d_vj < pratmo_ratio_limit & vmr_x_prat_dj <= vmr_limit_max(j));% & vmr_x_prat_dj >= vmr_limit_min(j); % the indexes of the columns of dj that aren't nans
                            %                         size(dj)
                            dj = dj(inonan); % remove nan values, and those that don't fit the criteria in the above line
                            %                         size(dj)
                            if ~isempty(dj)
                                mls_h_int_d_vj = reduce_tanstruct_by_rowindex(mls_h_int_d_vj, inonan); % reduce the data to the rows that aren't nans at that altitude
                                ratio_prat_h_int_d_vj = ratio_prat_h_int_d_vj(inonan); % reduce the ratios too. 1 x length(dj)
                                vmr_x_prat_dj = vmr_x_prat_dj(inonan); % reduce the vmrs * ratios too. 1 x length(dj)
                                [~,idjmin] = min(dj,[],2); % get the index of the smallest distance for altitude j.
                                %check if the value has been chosen before
                                chosen_rowcolumn_n = [j, find(mls.date_mjd == mls_h_int_d_vj.date_mjd(idjmin))]; % the row column of the chosen point in the main data (interpolated)
                                already_chosen = ismember(chosen_rowcolumn_n, chosen_rowcolumn, 'rows'); % this will be [1,1] if the value has already been chosen
                                if already_chosen == 1
                                    dj(idjmin) = nan; % change the value of dj to a nan and go back and choose another value
                                    redo = 1;
                                    disp('value already chosen, searching again')
                                    fprintf('n = %i. j = %i\n', n, j)
                                else
                                    redo = 0; % you have found a value that has not already been chosen so go ahead
                                    % fill in the values for the chosen index
                                    out.vmr(j,n)        = mls_h_int_d_vj.vmr(j,idjmin);
                                    out.vmr_error(j,n)  = mls_h_int_d_vj.vmr_error(j,idjmin);
                                    %                                 out.lat_tangent(j,n)  = mls3h_int_dj.lat_tangent(idjmin);
                                    %                                 out.lon_tangent(j,n)  = mls3h_int_dj.lon_tangent(idjmin);
                                    %                                 out.lat_mls(j,n)        = mls3h_int_dj.lat(j,idjmin);
                                    %                                 out.lon_mls(j,n)        = mls3h_int_dj.lon(j,idjmin);
                                    out.distance(j,n)   = dj(idjmin);
                                    out.date_mjd_mls(j,n) = mls_h_int_d_vj.date_mjd(idjmin);
%                                     out.time_diff(j,n) = out.date_mjd_mls(j,n) - out.date_mjd(n);
                                    % record the rows and altitudes of the chosen value to
                                    % make sure that it doesn't repeat
                                    k = k+1;
                                    chosen_rowcolumn(k,:) = chosen_rowcolumn_n;
                                    out.lst_ratio(j,n) = ratio_prat_h_int_d_vj(idjmin);
                                end
                            else
                                fprintf('no values meet the fractional and absolute limit criteria, n = %i. j = %i\n', n, j)
                            end
                        end
                    end
                else
                    fprintf('no data points lie within %f km that match the vortex criteria. n = %i.\n', distance_lim, n)
                    % the nth 'out' column will remain all NaNs in this case
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
%     if mod(n,1000) == 0 || n == 1
%         %         partialsave = sprintf('mlsstruct_acesample_%i_%i', starton, n);
%         mlsstruct_acesample_partialsave = out; %#ok<NASGU>
%         disp('partial save made without scaling applied')
%         save(savedest,'mlsstruct_acesample_partialsave');
%     else
%     end
end
out.vmr = out.vmr .* out.lst_ratio;
out.vmr_error = out.vmr_error .* out.lst_ratio;
out.ratio_applied = true;

smilesstruct_acesample = out;
fprintf('\nSaving sampled data with scaling applied to %s\n', savedest);
save(savedest,'smilesstruct_acesample')
disp('done')
%
toc
end

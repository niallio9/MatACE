function [ mlsstruct_acesample, chosen_rowcolumn] = sample_and_scale_mls_for_ace(mlsstruct_in, tanstruct_in, pratstruct, vmr_limit_in, do_day_night)
%A function to sample MLS data according to the time/lat/lon/alt of ACE
%measurements. The two closest times are sampled, mainly to have a choice
%of start-time when using chemical box model to scale the data at a later
%date.

% *INPUT*
%           mls_in: STRUCTURE - contains the gas specific MLS data. Can be
%           created with 'extract_mls_data.m'.
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
nocc = length(ace.occultation);
lalt_ace = length(ace.altitude_km(:,1));
prat = pratstruct; %#ok<NASGU>
pratmo_ratio_limit = Inf;
% pratmo_ratio_limit = pratmo_limit;
% vmr_limit = 2e-9;
vmr_limit = vmr_limit_in;
savedest = sprintf('MLS_v4p2_ClO_acesample_12h_1000km_vmrlimit%0.1d_20042010',vmr_limit) %#ok<NOPRT>

% lprat = length(prat.occultation);
% orbitprat = nan(2,lprat);
% orbitprat(1,1:lrat) = prat.occultation;
% orbitprat(2,1:lrat) = prat.sr1ss0;

%output structure
out.source_file = mls.source_file;
out.occultation = ace.occultation;
out.sr1ss0 = ace.sr1ss0;
out.beta_angle = nan(1,nocc);
out.date_mjd = nan(size(ace.vmr));
out.gas = mls.gas;
out.altitude_km = ace.altitude_km;
out.vmr = nan(size(ace.vmr));
out.vmr_error = nan(size(ace.vmr));
out.lat_tangent = nan(size(ace.vmr));
out.lon_tangent = nan(size(ace.vmr));
out.quality_flags = nan(size(ace.vmr));
out.pressure_hPa = ace.pressure_hPa;
out.lat = nan(size(ace.vmr));
out.lon = nan(size(ace.vmr));
out.distance = nan(size(ace.vmr));
out.time_diff = nan(size(ace.vmr));

chosen_rowcolumn = nan(lalt_ace.*nocc,2); % a matrix to hold the values of the chosen
% ratio_out = nan(lalt_ace, nocc);
out.lst_ratio = nan(lalt_ace, nocc);
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
    % subset mls data to 'time_lim' hours around the ace data
    Itime = find(abs(mls.date_mjd - ace_n.date_mjd) < time_lim);
    if ~isempty(Itime)
        %                 fprintf('\n%i MLS profiles within %0.2f days\n', length(Itime), time_lim)
        mls3h = reduce_tanstruct_by_rowindex(mls, Itime); % % get the mls data that is within time_lim of the ace time. I3h size data structure
        %% Add a new section to locate the biggest possible lat lon range and then subset the data. NJR - 07/2018
%         minlon_acen1 = min(ace_n.lon)
%         maxlon_acen1 = max(ace_n.lon)
%         minlat_acen1 = min(ace_n.lat)
%         maxlat_acen1 = max(ace_n.lat)
        
        minlon_acen = min(ace_n.lon) * (pi/180);
        maxlon_acen = max(ace_n.lon) * (pi/180);
        minlat_acen = min(ace_n.lat) * (pi/180);
        maxlat_acen = max(ace_n.lat) * (pi/180);
        % use the distance criteria to find a lat/lon range
        dlat = distance_lim ./ Re; % from l = r * dlat. the one for longitude is: l = r*cos(lat) * dlon
        dlon = distance_lim ./ (Re * cos(min(abs([minlat_acen, maxlat_acen])))); % use the lat that is closet to the equator
        %         dlon = 1.1*distance_lim ./ (Re * cos(minlat_acen)); % just use the minimum
        %         dlon_maxlat_acen = 1.1*distance_lim ./ (Re * cos(maxlat_acen)); % need two calculations for the lon range
        %         dlon_minlat_acen = 1.1*distance_lim ./ (Re * cos(minlat_acen));
        %         dlon = max([])
        latrange = [minlat_acen - dlat, maxlat_acen + dlat] * (180/pi);
        lonrange = [minlon_acen - dlon, maxlon_acen + dlon] * (180/pi);
        
        mls3h = subset_ace_by_lat(mls3h, latrange(1), latrange(2)); % subset to latitude range
        mls3h = subset_ace_by_lon(mls3h, lonrange(1), lonrange(2)); % subset to the longitude range
%         size(mls3h.occultation)
        
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
%                 size(mls3h_int_d.occultation)
                %% New section to alter which values to keep, depending on day or night LST
                if do_day_night == 1
                    [~,day_index] = subset_ace_by_lst_tangent(mls3h_int_d,6,18);
                    [~,night_index] = subset_ace_by_lst_tangent(mls3h_int_d,18,6);
                    press_highalt = find(mls3h_int_d.pressure_hPa(:,1) <= 30);
                    press_lowalt = find(mls3h_int_d.pressure_hPa(:,1) >= 30);
                    mls3h_int_d.vmr(press_highalt,day_index) = nan; % remove day values above 30hPa. values above 30hPa are to be from night time
                    mls3h_int_d.lat(press_highalt,day_index) = nan;
                    mls3h_int_d.lon(press_highalt,day_index) = nan;
                    mls3h_int_d.vmr(press_lowalt,night_index) = nan; % remove night values below 30hPa. values below 30hPa are to be from day time
                    mls3h_int_d.lat(press_lowalt,night_index) = nan;
                    mls3h_int_d.lon(press_lowalt,night_index) = nan;
                    %             zmls = p2z_waccm(mls3h_int_d.pressure_hPa.*100)/1000; % calculate the MLS altitude
                end
                %%
                lstace_n_rep = repmat(lstace_n, size(mls3h_int_d.vmr(1,:))); %#ok<NASGU> % duplicate the lst vector to be the same number of profiles as mls measurements
                %             whos
                %             size(mls3h_int_d.vmr(1,:))
                evalc('[~, ratio_prat3h_int_d] = scale_ace_with_pratmo_vmrs(mls3h_int_d, prat, lstace_n_rep);'); % get the pratmo ratios for the nth occultation and each of the mls profiles
                %             ratio_prat3h_int_d
                %             test = mls3h_int_d;
                %             return
                %need to recalculate the distances because the last reduction of the data above might
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
                    
                    mls3h_int_dj = mls3h_int_d; % initialise a structure that will be altered into a subset of the 'mls3h_int_d' in the loop below
                    vmr_x_prat_dj = mls3h_int_dj.vmr(j,:) .* ratio_prat3h_int_dj;
                    redo = 1; % to start the while loop
                    %                 if nansum(dj) == 0
                    %                     fprintf('no good values in time/distance criteria. n = %i. j = %i\n', n, j)
                    %                 end
                    %                 if nansum(ratio_prat3h_int_dj) == 0
                    %                     fprintf('no good pratmo values. n = %i. j = %i\n', n, j)
                    %                 end
                    %% main loop for finding a value at (n,j)
                    while nansum(dj) ~= 0 && nansum(ratio_prat3h_int_dj) ~= 0 && redo == 1
                        %                                                             disp('loop')
                        % % %                     inonan = find(~isnan(dj)); % the indexes of the columns of dj that aren't nans
                        inonan = find(~isnan(dj) & ~isnan(ratio_prat3h_int_dj) & ratio_prat3h_int_dj < pratmo_ratio_limit & abs(vmr_x_prat_dj) < vmr_limit); % the indexes of the columns of dj that aren't nans
                        dj = dj(inonan); % remove nan values
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
                                % % %                     elseif ratio_prat3h_int_dj(idjmin) > 100
                                % % %                         dj(idjmin) = nan; % change the value of dj to a nan and go back and choose another value
                                % % %                         redo = 1;
                                % % %                         disp('pratmo data is nan, searching again')
                                % % %                         fprintf('n = %i. j = %i\n', n, j)
                                % % %                     elseif ratio_prat3h_int_dj(idjmin) > 1
                                % % %                         dj(idjmin) = nan; % change the value of dj to a nan and go back and choose another value
                                % % %                         redo = 1;
                                % % %                         disp('pratmo ratio is > 1, searching again')
                                % % %                         fprintf('n = %i. j = %i\n', n, j)
                            else
                                redo = 0; % you have found a value that has not already been chosen so go ahead
                                % fill in the values for the chosen index
                                out.date_mjd(j,n)     = mls3h_int_dj.date_mjd(idjmin);
                                out.vmr(j,n)        = mls3h_int_dj.vmr(j,idjmin);
                                out.vmr_error(j,n)  = mls3h_int_dj.vmr_error(j,idjmin);
                                out.lat_tangent(j,n)  = mls3h_int_dj.lat_tangent(idjmin);
                                out.lon_tangent(j,n)  = mls3h_int_dj.lon_tangent(idjmin);
                                out.lat(j,n)        = mls3h_int_dj.lat(j,idjmin);
                                out.lon(j,n)        = mls3h_int_dj.lon(j,idjmin);
                                out.distance(j,n)   = dj(idjmin);
                                out.time_diff(j,n) = out.date_mjd(j,n) - ace_n.date_mjd;
                                
                                %                             out.distance(j,n)
                                %                             out.lat(j,n)
                                %                             ace_n.lat(j)
                                % record the rows and altitudes of the chosen value to
                                % make sure that it doesn't repeat
                                k = k+1;
                                chosen_rowcolumn(k,:) = chosen_rowcolumn_n;
                                %                     chosen_rowcolumn(1:k,:)
                                %                     ratio_out(j,n) = ratio_prat3h_int_dj(j, idjmin);
                                out.lst_ratio(j,n) = ratio_prat3h_int_dj(idjmin);
                            end
                        end
                    end
                    % % %                 %% fill in the values for the chosen index
                    % % %                 if ~nansum(dj) == 0 % if dj is not all nan values, either originally or because all values were already chosen
                    % % %                     out.date_mjd(j,n)     = mls3h_int_dj.date_mjd(idjmin);
                    % % %                     out.vmr(j,n)        = mls3h_int_dj.vmr(j,idjmin);
                    % % %                     out.vmr_error(j,n)  = mls3h_int_dj.vmr_error(j,idjmin);
                    % % %                     out.lat_tangent(j,n)  = mls3h_int_dj.lat_tangent(idjmin);
                    % % %                     out.lon_tangent(j,n)  = mls3h_int_dj.lon_tangent(idjmin);
                    % % %                     out.lat(j,n)        = mls3h_int_dj.lat(j,idjmin);
                    % % %                     out.lon(j,n)        = mls3h_int_dj.lon(j,idjmin);
                    % % %                     size(dj)
                    % % %                     out.distance(j,n)   = dj(idjmin);
                    % % %                     out.time_diff(j,n) = out.date_mjd(j,n) - ace_n.date_mjd;
                    % % %                     % record the rows and altitudes of the chosen value to
                    % % %                     % make sure that it doesn't repeat
                    % % %                     k = k+1;
                    % % %                     chosen_rowcolumn(k,:) = chosen_rowcolumn_n;
                    % % %                     %                     chosen_rowcolumn(1:k,:)
                    % % %                     %                     ratio_out(j,n) = ratio_prat3h_int_dj(j, idjmin);
                    % % %                     out.lst_ratio(j,n) = ratio_prat3h_int_dj(idjmin);
                    % % %                 else
                    % % %                     % the jth 'out' altitude for the nth occultation will remain a NaN in this case
                    % % %                 end
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

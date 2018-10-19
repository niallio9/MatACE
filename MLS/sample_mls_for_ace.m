function [mls_sample, chosen_rowcolumn] = sample_mls_for_ace(mlsstruct_in, tanstruct_in, do_day_night)
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
distance_lim = Inf; % the distance colocation criteria, in km.

if ~isfield(tanstruct_in,'lat')
    error('There is no GLC lat/lon information in the ACE structure. Stopping');
end
mls = convert_mls_to_ace_format(mlsstruct_in); % do this so that you can use the ace functions on the mls data
nocc = length(ace.occultation);
lalt_ace = length(ace.altitude_km(:,1));

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

chosen_rowcolumn = nan(lalt_ace.*nocc,2); % a matrix to hold the values of the chosen 
k = 0;

%% go through ace occultations
disp('going through ACE occultations to find coincident MLS points...')
for n = 1:nocc
    if ~rem(n,100)
        fprintf('past occultation %i of %i\n', n, nocc);
    end
    % subset mls data to 'time_lim' hours around the ace data
    Itime = find(abs(mls.date_mjd - ace.date_mjd(n)) < time_lim);
    if ~isempty(Itime)
%         fprintf('\n%i MLS profiles within %0.2f days\n', length(Itime), time_lim)
        mls3h = reduce_tanstruct_by_rowindex(mls, Itime); %#ok<NASGU> % get the mls data that is within time_lim of the ace time. I3h size data structure
        % interpolate all to the ace altitudes
        evalc('mls3h_int = interpolate_ace_to_pgrid(mls3h, ace.pressure_hPa(:,n));'); % interpolate the mls data to the pressure levels of the ace occultation. suppress command line output.
        % get the distance between the ace measurement and the mls
        % measurements, at each altitude, in km.
        d = latlon2distance(mls3h_int.lat,mls3h_int.lon,ace.lat(:,n), ace.lon(:,n), mls3h_int.altitude_km); % ace altitude x I3h array
        Id = find(d < distance_lim); % linear index here
        if ~isempty(Id)
%             fprintf('%i data points lie within %f km\n', length(Id), distance_lim)
            % subset to the data that is within 'distance_lim' of the ace measurement
            mls3h_int_d = reduce_tanstruct_data_by_index(mls3h_int, Id); % this will give an empty matrix if there are no points that fit the criteria.
            %% New scetion to alter which values to keep, depending on day or night LST
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
            %need to recalculate the distances because the last reduction of the data above might
            %have removed some measurements.
            d = latlon2distance(mls3h_int_d.lat, mls3h_int_d.lon,ace.lat(:,n), ace.lon(:,n), mls3h_int_d.altitude_km); % ace altitude x I3h_d array
            % this new d will have all values less than the distance limit.
            % choose the closest values at each altitude. Have to go row by
            % row because the 'min' function can't deal with nans properly.
            for j = 1:length(mls3h_int_d.altitude_km(:,1))
                % construct a new mls profile measurement from each of
                % the closest measurements
                dj = d(j,:); % the d values at the jth altitude
                redo = 1; % to start the while loop
                while nansum(dj) ~= 0 && redo == 1
                    inonan = find(~isnan(dj)); % the indexes of the columns of dj that aren't nans
                    dj = dj(inonan); % remove nan values
                    mls3h_int_dj = reduce_tanstruct_by_rowindex(mls3h_int_d, inonan); % reduce the data to the rows that aren't nans at that altitude
                    [~,idjmin] = min(dj,[],2); % get the index of the smallest distance for altitude j.
                    %check if the value has been chosen before
                    chosen_rowcolumn_n = [j, find(mls.date_mjd == mls3h_int_dj.date_mjd(idjmin))]; % the row column of the chosen point in the main data (interpolated)
                    already_chosen = ismember(chosen_rowcolumn_n, chosen_rowcolumn, 'rows'); % this will be [1,1] if the value has already been chosen
                    if already_chosen == 1 %sum(already_chosen) == 2
                        dj(idjmin) = nan; % change the value of dj to a nan and go back and choose another value
                        redo = 1;
                        disp('value already chosen, searching again')
                        fprintf('n = %i',n)
                    else
                        redo = 0; % you have found a value that has not already been chosen so go ahead
                    end
                end
                %% fill in the values for the chosen index
                if ~nansum(d(j,:)) == 0 % if dj is not all nan values, either originally or because all values were already chosen
                    out.date_mjd(j,n)     = mls3h_int_dj.date_mjd(idjmin);
                    out.vmr(j,n)        = mls3h_int_dj.vmr(j,idjmin);
                    out.vmr_error(j,n)  = mls3h_int_dj.vmr_error(j,idjmin);
                    out.lat_tangent(j,n)  = mls3h_int_dj.lat_tangent(idjmin);
                    out.lon_tangent(j,n)  = mls3h_int_dj.lon_tangent(idjmin);
                    out.lat(j,n)        = mls3h_int_dj.lat(j,idjmin);
                    out.lon(j,n)        = mls3h_int_dj.lon(j,idjmin);
                    % record the rows and altitudes of the chosen value to
                    % make sure that it doesn't repeat
                    k = k+1;
                    chosen_rowcolumn(k,:) = chosen_rowcolumn_n;
%                     chosen_rowcolumn(1:k,:)
                else
                    % the jth 'out' altitude for the nth occultation will remain a NaN in this case
                end
            end
        else
%             n
%             fprintf('no data points lie within %f km\n', distance_lim)
            % the nth 'out' column will remain all NaNs in this case
        end
    end
end
mls_sample = out;
disp('done')
%
toc
end

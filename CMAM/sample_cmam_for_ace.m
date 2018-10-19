function [cmam_sample] = sample_cmam_for_ace(filename_cmam, tanstruct_in)
%A function to sample cmam data according to the time/lat/lon/alt of ACE
%measurements. The two closest times are sampled, mainly to have a choice
%of start-time when using chemical box model to scale the data at a later
%date.

% *INPUT*    
%           filename_cmam: STRING - the name of the .nc file that holds the
%           cmam data.
%
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           cmam_sample: STRUCTURE - with the data that has been sampled
%           from CMAM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 01/18
%   NJR - 05/18 Now linearly interpolates in lat/lon/alt, and spline
%   interpolates in time
tic
%% Filter the ace data
ace = apply_ace_flags(tanstruct_in);
ace = filter_ace_pressure(ace);
ace = filter_ace_bad_lat(ace,10);
ace = filter_ace_bad_lon(ace,30);

%% Define some things
interptype = 'linear';
if ~isfield(tanstruct_in,'lat')
    error('There is no GLC lat/lon information in the ACE structure. Stopping');
end
cfile = filename_cmam;
cmam = read_cmam_ncdata(cfile); % a matlab structure with the data
nocc = length(ace.occultation);
lat_ace = ace.lat;
lon_ace = ace.lon; lon_ace(lon_ace < 0) = lon_ace(lon_ace < 0) + 360;
logp_ace = log(ace.pressure_hPa);
lat_cmam = cmam.lat;
lon_cmam = cmam.lon;
dlon_cmam = lon_cmam(2) - lon_cmam(1); % get the cmam lon step
logp_cmam = log(cmam.pressure_hPa);
%output structure
cout.source_file = cmam.source;
cout.occultation = ace.occultation;
cout.sr1ss0 = ace.sr1ss0;
cout.beta_angle = ace.beta_angle;
cout.date_mjd = nan(1,nocc);
% cout.date_mjd = nan(2, nocc);
cout.gas = cmam.gas;
cout.altitude_km = ace.altitude_km;
cout.vmr = nan(size(ace.vmr));
cout.vmr_error = zeros(size(ace.vmr)); % THIS MIGHT NEED TO BE CHANGED LATER
cout.lat_tangent = ace.lat_tangent;
cout.lon_tangent = ace.lon_tangent;
cout.quality_flags = nan(size(ace.quality_flags)); % THIS SHOULD BE CHANGED TO NANS IF YOU WANT TO GET ALL CMAM ALTITUDE LEVELS AND NOT WORRY ABOUT ACE FLAGS
cout.pressure_hPa = ace.pressure_hPa;
cout.lat = ace.lat;
cout.lon = ace.lon;
lmjd_cmam = length(cmam.date_mjd);
% cout.vmr1 = nan(size(ace.vmr));
% cout.vmr2 = nan(size(ace.vmr));

%% go through ace occultations
for n = 1:nocc
%     n
    if ~rem(n,1000)
        fprintf('\npast %i of %i\n', n, nocc)
    end
    lon_cmami = lon_cmam; % reset the lon grid each time
    [~,index] = sort(abs(ace.date_mjd(n)-cmam.date_mjd)); % sort the cmam times by distance from ace time
    ind_mjd1 = index(1); % get the index of the closest time
    if ind_mjd1 > 2 && ind_mjd1 < lmjd_cmam - 1
        ind_mjd5 = ind_mjd1 - 2 : ind_mjd1 + 2; % get two times either side of that
        mjd5_cmam = cmam.date_mjd(ind_mjd5);
    elseif ind_mjd1 <= 2
        ind_mjd5 = 1:5;
        mjd5_cmam = cmam.date_mjd(ind_mjd5);
    elseif ind_mjd1 >= lmjd_cmam - 1
        ind_mjd5 = lmjd_cmam - 4 : lmjd_cmam;
        mjd5_cmam = cmam.date_mjd(ind_mjd5);
    end
    %     ind_mjd2 = index(2);
    %     mjd_cmam12 = [cmam.date_mjd(ind_mjd1), cmam.date_mjd(ind_mjd2)]; % for interpolating in time
    %     vmrt1_cmam = squeeze(cmam.vmr(:,:,:,ind_mjd1)); % pull out the corresponding cmam vmrs: lon x lat x pres
    vmrt5_cmam = cmam.vmr(:,:,:,ind_mjd5); % pull out the corresponding cmam vmrs: lon x lat x pres x 5(time)
    %     vmrt2_cmam = squeeze(cmam.vmr(:,:,:,ind_mjd2));
    
    % Add a layer of longitude on either side of the vmr array to account
    % for interpolation near the limits of the lon grid. This creates the
    % illusion of a circular grid for linear interpolation.
    if max(lon_ace(:,n)) > lon_cmam(end)
        lon_cmami = cat(1,lon_cmam, lon_cmam(end) + dlon_cmam); % add another point to the end of the lon grid
        vmrt5_cmam = cat(1, vmrt5_cmam, vmrt5_cmam(1,:,:,:)); % add the first lon vmrs after the last
        %         vmrt2_cmam = cat(1, vmrt2_cmam, vmrt2_cmam(1,:,:));
    elseif min(lon_ace(:,n)) < lon_cmam(1)
        lon_cmami = cat(1, lon_cmam(1) - dlon_cmam, lon_cmam); % add a lon grid point before the first
        vmrt5_cmam = cat(1, vmrt5_cmam(end,:,:,:), vmrt5_cmam); % add the last lon vmrs before the first
        %         vmrt2_cmam = cat(1, vmrt2_cmam(end,:,:), vmrt2_cmam);
    end
    
    % interpolate these vmrs to the ACE lon/lat/alt
    %     vmrt1_int = interpn(lon_cmami, lat_cmam, logp_cmam, vmrt1_cmam, lon_ace(:,n), lat_ace(:,n), logp_ace(:,n), interptype, nan);
% % %     vmrt5_int = interpn(lon_cmami, lat_cmam, logp_cmam, mjd5_cmam, vmrt5_cmam, lon_ace(:,n), lat_ace(:,n), logp_ace(:,n), repmat(ace.date_mjd(n),150,1), interptype, nan);
    
    % linear interpolate in lat,lon,alt
    vmrt5_int = nan(length(logp_ace(:,n)),length(ind_mjd5));
    for i = 1:length(ind_mjd5)
        vmrt5_int(:,i) = interpn(lon_cmami, lat_cmam, logp_cmam, squeeze(vmrt5_cmam(:,:,:,i)), lon_ace(:,n), lat_ace(:,n), logp_ace(:,n), interptype, nan);
    end
    % spline interpolate in time
    vmrt1_int = nan(length(logp_ace(:,n)),1);
    for j = 1:length(logp_ace(:,n))
        if sum(~isnan(vmrt5_int(j,:))) < 4 % need 4 points for a spline interpolation
            vmrt1_int(j) = nan;
        else
            vmrt1_int(j) = interp1(mjd5_cmam, vmrt5_int(j,:)', ace.date_mjd(n), 'spline', nan);
        end
    end
    
    %     vmrt2_int = interpn(lon_cmami, lat_cmam, logp_cmam, vmrt2_cmam, lon_ace(:,n), lat_ace(:,n), logp_ace(:,n), interptype, nan);
    %     vmrt12_int(1,:,:) = vmrt1_int;
    %     vmrt12_int(2,:,:) = vmrt2_int;
    %     vmr_int = interp1(mjd_cmam12, vmrt12_int, ace.date_mjd(n));
    
    % fill in the output structure
    cout.vmr(:,n) = vmrt1_int;
    cout.date_mjd(1,n) = ace.date_mjd(n);
    %     cout.date_mjd(1,n) = cmam.date_mjd(ind_mjd1);
    %     cout.date_mjd(1,n) = cmam.date_mjd(ind_mjd1);
    %     cout.date_mjd(2,n) = cmam.date_mjd(ind_mjd2);
    %     cout.vmr1(:,n) = vmrt1_int;
    %     cout.vmr2(:,n) = vmrt2_int;
    %     cout.vmr(:,n) = vmr_int;
end
cmam_sample = cout;
%
toc
end

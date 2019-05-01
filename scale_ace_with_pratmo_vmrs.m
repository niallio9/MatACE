 function [ tanstruct_scaled, vmr_ratio ] = scale_ace_with_pratmo_vmrs( tanstruct_in, ratios_in, lst_in )
%A function to scale ACE VMR profiles to a given local solar time using
%VMRs calculated using a chemical box model (like PRATMO). PRATMO VMRs can
%be calculated with 'make_ace_gas_vmrs_with_pratmo.m'. When there are
%occultations present in the ace structure and not in the pratmo structure
%(and vice-versa), the vmr and error values are replaced with nans

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. The GLC or DMP data must also be
%           added to the tanstruct.
%
%           ratios_in: STRUCTURE - contains the vmr ratio data that has
%           been calculated using PRATMO.
%
%           lst_in: VECTOR - the local solar time(s) to which you want to
%           scale the ACE data.
%
% *OUTPUT*
%           tanstruct_scaled: STRUCTURE - with the same fields as the input
%           tanstruct, but with all of the volume mixing ratios scaled to
%           the input local solar time. The errors are also scaled using
%           the same factor. The time of the each measurement has been
%           edited to be the same day but at the local time of the input.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 03/18

%% define some things
%%USER DEFINED

%STANDARD
lstin = lst_in;
gas = tanstruct_in;
gasname = gas.gas;
rat = ratios_in;
ratname = rat.gas;
lzrat = length(rat.altitude_km(:,1));
lrat = length(rat.occultation);

fprintf('\nPROCESSING %s\n',gasname)
%check that the name of the gas for each structure is the same
if ~strcmp(ratname, gasname) % make sure that the name of the gas in the loaded file matches the one in the 
    warning('the name of the gas in the ACE structure (%s) does not match the name of the gas in the ratios structure (%s)', gasname, ratname);
end
%check that the length of the lst input matches the tanstruct input
if ~isscalar(lst_in) && isrow(lst_in) && length(lst_in) ~= length(gas.occultation)
    error('the length of the input LST vector should be the same as number of measurements to be scaled\n Maybe match the data first?')
end
if gas.altitude_km(1:lzrat,1) ~= rat.altitude_km(1:lzrat,1)
    error('the altitudes of the ACE and PRATMO data don''t match')
end

%% find where the occultations of the ace and pratmo data match
gas = apply_ace_flags(gas); % get rid of the -999 values, etc.
lgas = length(gas.occultation);
lzgas = length(gas.altitude_km(:,1));
% % % lstgas = get_ace_lst_tangent(gas);
lstgas = get_ace_lst(gas, 1); % lzgas x lgas
if isscalar(lstin)
%     lstin = lstin.*ones(1,lgas); % change to a vector of the same size as the ace input
    lstin = lstin.*ones(size(lstgas)); % change to a vector of the same size as the ace input
end

gasorbit = nan(2,lgas);
ratorbit = nan(2,lrat);
gasorbit(1,1:lgas) = gas.occultation;
gasorbit(2,1:lgas) = gas.sr1ss0;
ratorbit(1,1:lrat) = rat.occultation;
ratorbit(2,1:lrat) = rat.sr1ss0;
%Check which ones match
% [~,ygas,yrat] = intersect(gasorbit',ratorbit','rows'); % the indices of where the orbits/occultations match
[ygas, yrat] = ismember(gasorbit',ratorbit','rows'); % use this instead of above line so that repeated occultations can be used
[ygas] = find(ygas);
[~,~,yrat] = find(yrat);
% lygas = length(ygas);
% x = ratorbit;
% y = gasorbit;
% return

vmr_lst_ace = nan(lzgas,lgas);
vmr_lst_in = nan(lzgas,lgas);
% vmr_scaled = nan(lzgas,lgas);
% vmr_error_scaled = nan(lzgas,lgas);

%% set up new tanstruct
out.source_file = gas.source_file;
out.occultation = nan(1,lgas);
out.sr1ss0 = nan(1,lgas);
out.beta_angle = nan(1,lgas);
out.date_mjd = nan(lzgas,lgas);
% out.date_mjd = rat.date_mjd(yrat); % change the time to the times of the ratios. This will be the same day but at the 'LST_in'
out.gas = gasname;
out.altitude_km = gas.altitude_km;
out.vmr = nan(lzgas,lgas);
out.vmr_error = nan(lzgas,lgas);
out.lat_tangent = nan(1,lgas);
out.lon_tangent = nan(1,lgas);
out.quality_flags = nan(lzgas,lgas);
out.pressure_hPa = nan(lzgas,lgas);
if isfield(gas,'lon')
    out.lon = nan(lzgas,lgas);
    out.lat = nan(lzgas,lgas);
end

%% get the PRATMO profiles at lst_ace and lst_in and find the ratio by VMR@lst_in/VMR2lst_ace
disp('scaling data...')
for i = 1:length(ygas)
%         i
    inonan = ~isnan(rat.lst(:,yrat(i))); % get indices of where there are not nans in the pratmo LSTs
    if sum(inonan) > 10 % only do if you have more than 10 profiles in pratmo for that day
        for j = 1:lzrat
            ratlsti = rat.lst(inonan,yrat(i)); % inonan x 1
            ratvmri = rat.vmr(inonan,j,yrat(i)); % inonan x 1 x 1
            [ratlsti, Iratlsti] = unique(ratlsti); % remove repeating lst values for this profile
            ratvmri = ratvmri(Iratlsti); % remove the corresponding vmrs too
            vmr_lst_ace(j,ygas(i)) = interp1(ratlsti, ratvmri, lstgas(j, ygas(i))); % lzgas x lgas
            vmr_lst_in(j,ygas(i)) = interp1(ratlsti, ratvmri, lstin(j, ygas(i))); % lzgas x lgas
        end
    end
end 

%% scale and fill in the output data
if ~isempty(ygas)
    vmr_ratio = vmr_lst_in ./ vmr_lst_ace; % lzgas x lgas. has nans at occultations not in ygas
    vmr_ace = gas.vmr; % lzgas x lgas
    
    ifound = find(vmr_ace < 0 & vmr_ratio > 1); % vmr < 0, pratmo > 1
    if ~isempty(ifound)
        %                         disp('making some negative vmrs positive for >1 scaling.')
        %                         fprintf('n = %i. j = %i\n', n, j)
        vmr_ace(ifound) = vmr_ace(ifound) * -1; % change the vmr from negative to positive here. this is the same as x + 2*abs(x)
        vmr_ace(ifound) = ( (vmr_ratio(ifound) .* vmr_ace(ifound)) - 2*vmr_ace(ifound) ) ./ vmr_ratio(ifound);
        % this line above gives:
        % ((pratmo*vmr) - 2*vmr) / pratmo
        % which, when multiplied by pratmo later on,
        % will give (pratmo*vmr) - 2*vmr)
    end
    
    
    vmr_scaled = vmr_ace .* vmr_ratio; % lzgas x lgas. nans at occultations not in ygas
    vmr_error_scaled = gas.vmr_error .* vmr_ratio; % lzgas x lgas. 
    % vmr_scaled(1:lzrat,:) = gas.vmr(1:lzrat,ygas) .* vmr_ratio;
    % vmr_error_scaled(1:lzrat,:) = gas.vmr_error(1:lzrat,ygas) .* vmr_ratio;
else % otherwise give out nans the same length as the number of measurements
    warning('there are no matching occultations')
    vmr_ratio = nan(lzgas,lgas);
    vmr_scaled = vmr_ratio; % nans
    vmr_error_scaled = vmr_ratio; % nans
end
%%output the new tanstruct
out.source_file = gas.source_file;
out.occultation = gas.occultation;
out.sr1ss0 = gas.sr1ss0;
out.beta_angle = gas.beta_angle;
% whos
% size(gas.date_mjd)
if isvector(gas.date_mjd) && ismatrix(lstgas)
    gas.date_mjd = repmat(gas.date_mjd, length(lstgas(:,1)), 1); % to run on stupid fucking deluge
end
% length(ygas)
out.date_mjd = gas.date_mjd + (lstin - lstgas)./24; % change the time to the same day but at the lst_in. Shift the original time by the number of hours between the original and new LST
% out.date_mjd = rat.date_mjd(yrat); % change the time to the times of the ratios. This will be the same day but at the 'LST_in'
out.gas = gasname;
out.altitude_km = gas.altitude_km;
out.vmr = vmr_scaled;
out.vmr_error = vmr_error_scaled;
out.lat_tangent = gas.lat_tangent;
out.lon_tangent = gas.lon_tangent;
% out.quality_flags = gas.quality_flags(:,ygas);
out.pressure_hPa = gas.pressure_hPa;
% if length(gas.pressure_hPa(1,:)) > 1
%     out.pressure_hPa = gas.pressure_hPa(:,ygas);
% else
%     out.pressure_hPa = gas.pressure_hPa;
% end
if isfield(gas,'lon')
    out.lon = gas.lon;
    out.lat = gas.lat;
end
out.lst_ratio = vmr_ratio;
out.ratio_applied = true;
tanstruct_scaled = out;

disp('done')
%
end


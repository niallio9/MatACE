function [ tanstruct_scaled ] = scale_ace_with_pratmo( tanstruct_in, tanstruct_o3_in, tanstruct_T_in, lst_in )
%A function to scale ACE VMR profiles to a given local solar time using a
%chemical box model.
% *****************NOTE, THIS FUNCTION IS OBSOLETE AND TAKES A VERY LONG
% TIME TO RUN. IT IS BETTER TO USE SOMETHING LIKE
% 'make_ace_nitrogen_ratios_with_pratmo.m' or
% 'make_ace_gas_vmrs_with_pratmo', WHICH CALCULATE THE RATIOS FOR ALL
% POSSIBLE ACE MEASUREMENTS AND FOR MULTIPLE GASES AT ONCE. THE RATIO
% INFORMATION IS THEN SAVED FOR LATER USE***************** NJR - 04/18

% *INPUT*    
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. The GLC data must also be added to
%           the tanstruct so that it has the latitude information.
%
%           tanstruct_O3_in: STRUCTURE - contains the O3 specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           tanstruct_T_in: STRUCTURE - contains the T specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           lst_in: FLOAT - the local solar time of the data
%
% *OUTPUT*
%           tanstruct_scaled: STRUCTURE - with the same fields as the input
%           tanstruct, but with all of the volume mixing ratios scaled to
%           the input local solar time. The errors are also scaled using
%           the same factor. The time of the each measurement has been
%           edited to be the same day but at the local time of the input.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 03/18

% %% make sure the ACE GLC info is included in the tanstruct
% %this can be added to the structure using 'merge_ace_glc.m'
% if ~isfield(tanstruct_o3_in,'lat')
%     error('There is no GLC lat/lon information in the ACE structure. Stopping');
% end

%% define some things
acegas = tanstruct_in;
gasname = acegas.gas;
if length(gasname) > 7
    if strcmp(gasname(end-5:end), 'sap_am') || strcmp(gasname(end-5:end), 'sap_pm')
        gasname_short = gasname(1:end-7);
        if lst_in < 12
            gasname_out = sprintf('%s_sap_s%02.0fam', gasname_short,lst_in);
        else
            gasname_out = sprintf('%s_sap_s%02.0fam', gasname_short,lst_in);
        end
    end
elseif length(gasname) > 3
    if strcmp(gasname(end-1:end), 'am') || strcmp(gasname(end-1:end), 'pm')
        gasname_short = gasname(1:end-3);
        if lst_in < 12
            gasname_out = sprintf('%s_s%02.0fam', gasname_short,lst_in);
        else
            gasname_out = sprintf('%s_s%02.0fam', gasname_short,lst_in);
        end;
    end
else
    if lst_in < 12
        gasname_out = sprintf('%s_s%02.0fam', gasname,lst_in);
    else
        gasname_out = sprintf('%s_s%02.0fam', gasname,lst_in);
    end
end
% gasname = 'NO2';
ifig = 0; %#ok<NASGU>
%% Filter the ace data (it's also done in subcode but it doesn't matter)
acegas = apply_ace_flags(acegas);
acegas = filter_ace_pressure(acegas);
% acegas = filter_ace_bad_lat(acegas,10);
% acegas = filter_ace_bad_lon(acegas,30);

aceo3 = apply_ace_flags(tanstruct_o3_in);
aceT = apply_ace_flags(tanstruct_T_in);
% don't need to do these because the
% two ace datasets are matched later on -NJR
% aceT = filter_ace_pressure(aceT); 
% aceT = filter_ace_bad_lat(aceT,10);
% aceT = filter_ace_bad_lon(aceT,30);

%% match the ace data
[aceo3, aceT] = match_ace_data(aceo3, aceT);
[aceo3, acegas] = match_ace_data(aceo3, acegas);
[aceT, acegas] = match_ace_data(aceT, acegas); % now all of the data should be matching 
%% get the ACE variables
nace = length(aceo3.occultation);
zlowi = aceo3.altitude_km < 75; % get data below 75km. for the box model, i think
% zlowi(1:15) = 0;
zlow = aceo3.altitude_km(zlowi);
lzlow = length(zlow);
o3 = aceo3.vmr(zlowi,:);
T = aceT.vmr(zlowi,:);
P = aceo3.pressure_hPa(zlowi,:);
kb = 1.38064852e-23; % Boltzmann's constant
Nair = (P./(T * kb))*1e-6; % this is the air density in cm^3
o3 = o3.*Nair;
lat = aceo3.lat_tangent;
lon = aceo3.lon_tangent;
lon(lon<0) = lon(lon<0) + 360;
jday = mjd2doy(aceo3.date_mjd);
lst_ace = get_ace_lst_tangent(aceo3);
lst_out = lst_in; % the time to which you want to scale the data
NOy = zeros(lzlow,1);
N2O = zeros(lzlow,1);
Bry = zeros(lzlow,1);
% outputs
vmr_scaled = nan(size(acegas.vmr));
vmr_error_scaled = nan(size(acegas.vmr));

%% load the albedo data
data = load('MERIS_412nm_albedo_avg.mat','data'); data = data.data;
alat = data.lat; alon = data.lon; alb = data.alb; clear data;
alon(alon<0) = alon(alon<0) + 360;

%% loop through the ace data and run the box model for each measurement
for j = 1:10%nace
    tic
    if mod(j,200) == 0
       fprintf('\noccultation %i of %i\n',j,nace) 
    end
    % set up the box model input
    init.lat = lat(j);
    init.jday = jday(j);
    latj = lat(j);
    lonj = lon(j);
    [~, ltI] = min(abs(alat-latj));
    [~, lnI] = min(abs(alon-lonj));
    a = alb(ltI, lnI);
    init.albedo = a;
    init.ndmr = 0;  % 0=save as vmr; 1=number density
    init.atmos = [zlow T(:,j) Nair(:,j) o3(:,j) NOy N2O Bry];
%     init.atmos = [];
    init.aero_scalefactor = 1;
    init.diurnalO3 = 0;
    init.z = zlow';   % model will provide output for radicals from about 10-58 km (zeros outside that)

    % run with O3 held constant
    [box0] = run_boxmodel_njr(init);
    
    %find the profile values at the LSTs of ACE and get the
    %scaling fraction
    lst_pratmo = box0.LST; % this is 34 times from 12noon to 12noon: 34x1 array
    [lst_pratmo, Ilst_pratmo] = sort(lst_pratmo(1:end-1)); % remove the second noon value and sort the values from 00:00 onwards. 33x1.
    lst_pratmo = [lst_pratmo ; 24]; %#ok<AGROW> % add a value for the second midnight, i.e., 24:00. 34x1
    Ilst_pratmo = [Ilst_pratmo ; Ilst_pratmo(1)]; %#ok<AGROW> % add the index of the second midnight value, which is the same as the first index (for 00:00). 34x1
    
    gas = eval(sprintf('box0.%s(1:end-1,:);', gasname_short)); % get the gas vmr data from the pratmo output
    gas = gas(Ilst_pratmo,:); % sort the gas data from 00:00 to 24:00. 34x75
    gas_lstace = interp1(lst_pratmo, gas, lst_ace(j));
    gas_lst_out = interp1(lst_pratmo, gas, lst_out);
    gas_rat = gas_lst_out' ./ gas_lstace'; % the ratio of the vmr at the new time and the vmr at the original ace time

    
    %% scale the ace gas data
    vmr_scaled(1:lzlow,j) = acegas.vmr(1:lzlow,j) .* gas_rat;
    vmr_error_scaled(1:lzlow,j) = acegas.vmr_error(1:lzlow,j) .* gas_rat;
%     %%some plotting
%     ifig = ifig + 1;
%     zstart = 1;
%     figure(ifig), plot(vmr_scaled(zstart:length(zlow),j),zlow(zstart:end),'b', acegas.vmr(:,j),acegas.altitude_km,'k')
%     title(sprintf(' ACE %s scaled', gasname),'interpreter','none')
%     xlabel(sprintf('%s scaled VMR',gasname),'interpreter','none')
%     ylabel('altitude [km]')
%     legend(sprintf('ACE scaled to %0.2f LST',lst_out), sprintf('ACE @ %0.2f LST',lst_ace(j)) )
%
%     figure(ifig+10), plot(gas_rat(zstart:end),zlow(zstart:end))
%     title(sprintf('ratio of %s at two LSTs', gasname),'interpreter','none')
%     xlabel(sprintf('%s ratio [VMR/VMR]',gasname),'interpreter','none')
%     ylabel('altitude [km]')
toc
end
%%output the new tanstruct
out = acegas;
out.date_mjd = aceo3.date_mjd + (lst_in - lst_ace)./24; % change the time to the same day but at 10am. Shift the original time by the number of hours between the origianl and new LST
out.gas = gasname_out;
out.vmr = vmr_scaled;
out.vmr_error = vmr_error_scaled;
tanstruct_scaled = out;
%
end


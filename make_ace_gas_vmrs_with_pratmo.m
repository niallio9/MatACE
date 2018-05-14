function [ ] = make_ace_gas_vmrs_with_pratmo( tanstruct_o3_in, tanstruct_T_in, varargin )
%A function to calculate the the VMRs of trace species at multiple local
%solar times (LSTs) throughout the days of ACE measurements, at the
%measurement locations. The values can be used to create ratios that can be
%used later to scale the ace measurements.

% *INPUT*    
%           tanstruct_O3_in: STRUCTURE - contains the O3 specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           tanstruct_T_in: STRUCTURE - contains the T specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           varargin: STRING - the name of the gas for which you want to
%           calculate VMR values with PRATMO.
%           The input may be a single gas (e.g., 'O3'), or multiple
%           gases (e.g., 'O3', 'ClO').
%
% *OUTPUT*
%           tanstruct_scaled: STRUCTURE - the VMR information as a function
%           of LST is saved for each input species in the current
%           directory.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 04/18

%% define some things

%USER DEFINED


%STANDARD
gasnames = varargin; % array of cells with gas names
lgases = length(gasnames); % number of gases to read
aceo3 = tanstruct_o3_in;
aceT = tanstruct_T_in;
aceo3 = apply_ace_flags(aceo3);
aceT = apply_ace_flags(aceT);

%% match the ace data
[aceo3, aceT] = match_ace_data(aceo3, aceT);

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
% lst_ace = get_ace_lst_tangent(aceo3);
% lst_out = lst_in; % the time to which you want to scale the data
NOy = zeros(lzlow,1);
N2O = zeros(lzlow,1);
Bry = zeros(lzlow,1);
% outputs
for i = 1:lgases
    setgas = sprintf('gasout_%s = nan(34,lzlow,nace);', gasnames{i});
    eval(setgas);
end
lst_pratmo_out = nan(34,nace);

%% load the albedo data
data = load('MERIS_412nm_albedo_avg.mat','data'); data = data.data;
alat = data.lat; alon = data.lon; alb = data.alb; clear data;
alon(alon<0) = alon(alon<0) + 360;

%% loop through the ace data and run the box model for each measurement
fprintf('\nRunning PRATMO for %i occultations...',nace)
tic
for j = 1:nace
    if mod(j,10) == 0
        fprintf('\npast %i of %i\n', j, nace)
        %        fprintf('\noccultation %i of %i milestone\n',j,nace)
        %        toc % check how long that took
        %        tic % start a new timer
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
    try % sometimes PRATMO won't converge. Skip this part if it doesn't
        lst_pratmo = box0.LST; % this changes in length sometimes. I don't know why. usually this is 34 times from 12noon to 12noon: 34x1 array, usually
        [lst_pratmo, Ilst_pratmo] = sort(lst_pratmo(1:end-1)); % remove the second noon value and sort the values from 00:00 onwards. 33x1, usually
        lst_pratmo = [lst_pratmo ; 24]; %#ok<AGROW> % add a value for the second midnight, i.e., 24:00. 34x1, usually
        Ilst_pratmo = [Ilst_pratmo ; Ilst_pratmo(1)]; %#ok<NASGU,AGROW> % add the index of the second midnight value, which is the same as the first index (for 00:00). 34x1, usually
        lpratmo_lst = length(lst_pratmo);
        lst_pratmo_out(1:lpratmo_lst,j) = lst_pratmo;
        
%         no = box0.NO(1:end-1,:) % 33x75, usually. first dimension changes sometimes. I don't know why.
%         whos
        for n = 1:lgases
            getgas_n = sprintf('gas_%s = box0.%s(1:end-1,:);', gasnames{n}, gasnames{n}); % 33x75 (lpratmo_lst-1 x 75)
            eval(getgas_n);
            sortgas_n = sprintf('gas_%s = gas_%s(Ilst_pratmo,:);', gasnames{n}, gasnames{n}); % 34X75
            eval(sortgas_n);
            fillgas_n = sprintf('gasout_%s(1:lpratmo_lst,:,j) = gas_%s;', gasnames{n}, gasnames{n}); % 34x75 into 34x75xnace.
            eval(fillgas_n);
%             gasout(1:lpratmo_lst-1,:,1) = box0.NO(1:end-1,:); % 33x75 (lpratmo_lst-1 x 75) into 34x75x6.
        end
    catch
    end
end
disp('done')
clear aceT o3 Nair P T alb gases jday lat lon n

%% save the output data
for n = 1:lgases
    %%output the new pratstruct
    out.occultation = aceo3.occultation;
    out.sr1ss0 = aceo3.sr1ss0;
%     out.beta_angle = aceo3.beta_angle;
    out.date_mjd = aceo3.date_mjd;
    out.gas = gasnames{n};
    out.altitude_km = zlow;
    setvmr_n = sprintf('out.vmr = gasout_%s;', gasnames{n}); % 34x75xnace
    eval(setvmr_n);
    out.lat_tangent = aceo3.lat_tangent;
    out.lon_tangent = aceo3.lon_tangent;
    out.lst = lst_pratmo_out; % same for all gases
    % out.lon = aceo3.lon;
    % out.lat = aceo3.lat;
    pratstruct = out; %#ok<NASGU>
    
    % save for now in case matlab closes or some shit
    savedest = sprintf('ACE_v3p6_pratmo_%s_all_LST',gasnames{n});
    fprintf('saving %s data to %s\n', gasnames{n}, savedest);
    save(savedest,'pratstruct')
    
end
disp('All done :)')
%
end


function [ ratios ] = make_ace_nitrogen_ratios_with_pratmo( tanstruct_o3_in, tanstruct_T_in, lst_input )
%A function to calculate the ratio of the VMRs of nitrogen-containing
%species at the input local solar time (LST) and at the ACE measurement
%times. The ratios can be used later to scale the ace measurements of the
%nitrogen species wit, e.g., 'scale_ace_with_nitrogen_ratios'.

% *INPUT*    
%           tanstruct_O3_in: STRUCTURE - contains the O3 specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           tanstruct_T_in: STRUCTURE - contains the T specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           lst_input: FLOAT - the local solar time of the data
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
% %this can be added to the structure using 'merge_ace_glc.m'. It should be
% %there because the input data here should have been subset by LST.
% %Subsetting can be done with 'subset_ace_by_lst'.
% if ~isfield(tanstruct_o3_in,'lat')
%     error('There is no GLC lat/lon information in the ACE structure. Stopping');
% end

%% define some things

%USER DEFINED
lst_in = lst_input; % changed this so that it is an input to the function

%STANDARD
if lst_in < 12
    fprintf('\nchosen LST < 12. Subsetting data to AM times\n')
    [aceo3,~] = split_ace_by_lst_tangent(tanstruct_o3_in);
    [aceT,~] = split_ace_by_lst_tangent(tanstruct_T_in);
else
    fprintf('\nChosen LST > 12. Subsetting data to PM times\n')
    [~, aceo3] = split_ace_by_lst_tangent(tanstruct_o3_in);
    [~, aceT] = split_ace_by_lst_tangent(tanstruct_T_in);
end
aceo3 = apply_ace_flags(aceo3);
aceT = apply_ace_flags(aceT);
gasnames = {'NO','NO2','ClONO2','N2O5','HNO3'};
lgases = length(gasnames);

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
Nair = ( 100*P ./ ( (T + 273.15) * kb) )*1e-6; % this is the air density in cm^3
o3 = o3.*Nair;
lat = aceo3.lat_tangent;
lon = aceo3.lon_tangent;
lon(lon<0) = lon(lon<0) + 360;
jday = mjd2doy(aceo3.date_mjd);
lst_ace = get_ace_lst_tangent(aceo3);
% lst_out = lst_in; % the time to which you want to scale the data
NOy = zeros(lzlow,1);
N2O = zeros(lzlow,1);
Bry = zeros(lzlow,1);
% outputs
gases = nan(34,lzlow,lgases); % pratmo outputs profiles at 34 LSTs
gases_rat = nan(lzlow,lgases,nace);

%% load the albedo data
data = load('MERIS_412nm_albedo_avg.mat','data'); data = data.data;
alat = data.lat; alon = data.lon; alb = data.alb; clear data;
alon(alon<0) = alon(alon<0) + 360;

%% loop through the ace data and run the box model for each measurement
fprintf('\nRunning PRATMO for %i occultations...',nace)
tic
for j = 1:nace
    if mod(j,10) == 0
        fprintf('\n%i\n',j)
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
        Ilst_pratmo = [Ilst_pratmo ; Ilst_pratmo(1)]; %#ok<AGROW> % add the index of the second midnight value, which is the same as the first index (for 00:00). 34x1, usually
        lpratmo_lst = length(lst_pratmo);
        
        %     no = box0.NO(1:end-1,:) % 33x75, usually. first dimension changes sometimes. I don't know why.
        %     no2 = box0.NO2(1:end-1,:);
        %     clono2 = box0.ClONO2(1:end-1,:);
        %     n2o5 = box0.N2O5(1:end-1,:);
        %     hno3 = box0.HNO3(1:end-1,:);
        %     whos
        gases(1:lpratmo_lst-1,:,1) = box0.NO(1:end-1,:); % 33x75 (lpratmo_lst-1 x 75) into 34x75x6.
        gases(1:lpratmo_lst-1,:,2) = box0.NO2(1:end-1,:);
        gases(1:lpratmo_lst-1,:,3) = box0.ClONO2(1:end-1,:);
        gases(1:lpratmo_lst-1,:,4) = box0.N2O5(1:end-1,:);
        gases(1:lpratmo_lst-1,:,5) = box0.HNO3(1:end-1,:);
        
        gases = gases(Ilst_pratmo,:,:); % sort the gas data from 00:00 to 24:00. 34x75x6
        gases_lstace = squeeze(interp1(lst_pratmo, gases, lst_ace(j))); % 75x6
        gases_lst_out = squeeze(interp1(lst_pratmo, gases, lst_in)); % 75x6
        gases_rat(:,:,j) = gases_lst_out ./ gases_lstace; % the ratio of the vmr at the new time and the vmr at the original ace time 75x6xnace
    catch
    end
end

clear aceT o3 Nair P T alb gases jday lat lon
%%output the new tanstruct
out.occultation = aceo3.occultation;
out.sr1ss0 = aceo3.sr1ss0;
out.beta_angle = aceo3.beta_angle;
out.date_mjd = aceo3.date_mjd + (lst_in - lst_ace)./24; % change the time to the same day but at 10am. Shift the original time by the number of hours between the origianl and new LST
out.gas = gasnames;
out.altitude_km = zlow;
out.vmr_ratio = permute(gases_rat,[2,1,3]); % change to a 5x75xnace
out.lat_tangent = aceo3.lat_tangent;
out.lon_tangent = aceo3.lon_tangent;
out.LST = lst_in;
% out.lon = aceo3.lon;
% out.lat = aceo3.lat;
ratios = out;
disp('All done :)')

% save for now in case matlab closes or some shit
savedest = sprintf('ACE_v3p6_pratmo_nitrogen_ratios_%2.0fLST_2018',lst_in);
fprintf('saving nitrogen data to %s\n', savedest);
save(savedest,'ratios')
%
end


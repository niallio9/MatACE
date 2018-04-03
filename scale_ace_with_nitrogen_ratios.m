function [ ] = scale_ace_with_nitrogen_ratios( ratios_in )
%A function to scale ACE VMR profiles to a given local solar time using VMR
%ratios calculated using a chemical box model (like PRATMO). Ratios can be
%calculated using 'make_ace_nitrogen_ratios_with_pratmo'.

% *INPUT*    
%           ratios_in: STRUCTURE - contains the vmr ratio data that has
%           been calculated using PRATMO.
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
home_linux = '/home/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\Dropbox\climatology\nryan\'; %#ok<NASGU>
home_deluge = '/net/deluge/pb_1/users/nryan/'; %#ok<NASGU>

% matdirectory = '/Volumes/Seagate Backup Plus Drive/ACE/matdata/';
% matdirectory = 'F:\ACE\matdata\';
% matdirectory = strcat(home_mac,'matdata/');
% matdirectory = '/Users/niall/ACE/matdata';
matdirectory = 'C:\Users\ryann\ACE\matdata';
% matdirectory = strcat(home_deluge,'ACE/','matdata/');
if ~isdir(matdirectory)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',matdirectory)
    error('The directory containing the .mat data couldn''t be found')
end
rat = ratios_in;
lst_in  = rat.LST;
% CHANGE THE GASES IN HERE IF YOU WANT TO SCALE DIFFERENT TYPES OF THE FIVE
% NITROGEN CLIMATOLOGIES 
if lst_in < 12
    gasnames = {'NO_sap_am', 'NO2_sap_am', 'ClONO2_sap_am', 'N2O5_sap_am', 'HNO3_sap_am' }; % an array with the gas names
else
    gasnames = {'NO_sap_pm', 'NO2_sap_pm', 'ClONO2_sap_pm', 'N2O5_sap_pm', 'HNO3_sap_pm' }; % an array with the gas names;
end
%STANDARD
lzrat = length(rat.altitude_km(:,1));
lrat = length(rat.occultation);
lst_out_name = lst_in - 12; % just for naming the output files


%% Go throught the gases, load the files, match the data, scale the data, save the scaled data
filein_pre = 'ACE_v3p6_';
filein_post = '.mat';

for i = 1:length(gasnames)
    fprintf('\nPROCESSING %s\n',gasnames{i})
    gasname = gasnames{i};
    filename = strcat(filein_pre,gasname);
    filein_i = fullfile(matdirectory, filename);
    load(filein_i); gas = tanstruct; clear tanstruct;
    if ~strcmp(gasname, gas.gas) % make sure that the name of the gas in the loaded file is the same as in the filename
        error('the name of the gas in %s does not match %s', filein_i, gasname);
    end
    %check that the name of the gas for each structure is the same
    ratname = rat.gas{i};
    if ~strcmp(ratname, gasname(1:length(ratname))) % make sure that the name of the gas in the loaded file matches the correct dimesion of the ratios array
        error('the name of the gas in %s (%s) does not match the name of the dimension of the ratios array (%s)', filein_i, gasname, ratname);
    end
    gas = apply_ace_flags(gas); % get rid of the -999 values, etc.
    lgas = length(gas.occultation);
    lzgas = length(gas.altitude_km(:,1));
    % get the output name
    if length(gasname) > 7 && strcmp(gasname(end-5:end-3), 'sap')
        if strcmp(gasname(end-5:end), 'sap_am') || strcmp(gasname(end-5:end), 'sap_pm')
            gasname_short = gasname(1:end-7);
            if lst_in < 12
                gasname_out = sprintf('%s_sap_s%02.0fam', gasname_short,lst_in);
            else
                gasname_out = sprintf('%s_sap_s%02.0fpm', gasname_short,lst_out_name);
            end
        end
    elseif length(gasname) > 3
        if strcmp(gasname(end-1:end), 'am') || strcmp(gasname(end-1:end), 'pm')
            gasname_short = gasname(1:end-3);
            if lst_in < 12
                gasname_out = sprintf('%s_s%02.0fam', gasname_short,lst_in);
            else
                gasname_out = sprintf('%s_s%02.0fam', gasname_short,lst_out_name);
            end;
        end
    else
        if lst_in < 12
            gasname_out = sprintf('%s_s%02.0fam', gasname,lst_in);
        else
            gasname_out = sprintf('%s_s%02.0fpm', gasname,lst_out_name);
        end
    end
    %%match the ace data
    gasorbit = nan(2,lgas);
    ratorbit = nan(2,lrat);
    gasorbit(1,1:lgas) = gas.occultation;
    gasorbit(2,1:lgas) = gas.sr1ss0;
    ratorbit(1,1:lrat) = rat.occultation;
    ratorbit(2,1:lrat) = rat.sr1ss0;
    %Check which ones match
    [~,ygas,yrat] = intersect(gasorbit',ratorbit','rows'); % the indices of where the orbits/occultations match
    vmr_scaled = nan(lzgas,length(ygas));
    vmr_error_scaled = nan(lzgas,length(ygas));
    %%scale the gas data with the ratio data. the ratios are
    %%VMR@LST_in / VMR@LST_ace
%     size(gas.vmr)
%     lzrat
%     whos
%     gasvmr = gas.vmr(1:lzrat,ygas); % for testing
%     ratvmr = squeeze(rat.vmr_ratio(i,:,yrat)); % for testing
%     whos
    vmr_scaled(1:lzrat,:) = gas.vmr(1:lzrat,ygas) .* squeeze(rat.vmr_ratio(i,:,yrat));
    vmr_error_scaled(1:lzrat,:) = gas.vmr_error(1:lzrat,ygas) .* squeeze(rat.vmr_ratio(i,:,yrat));
    
    %%output the new tanstruct
    out.source_file = gas.source_file;
    out.occultation = gas.occultation(ygas);
    out.sr1ss0 = gas.sr1ss0(ygas);
    out.beta_angle = gas.beta_angle(ygas);
    out.date_mjd = rat.date_mjd(yrat); % change the time to the times of the ratios. This will be the same day but at the 'LST_in'
    out.gas = gasname_out;
    out.altitude_km = gas.altitude_km;
    out.vmr = vmr_scaled;
    out.vmr_error = vmr_error_scaled;
    out.lat_tangent = gas.lat_tangent(ygas);
    out.lon_tangent = gas.lon_tangent(ygas);
    out.quality_flags = gas.quality_flags(:,ygas);
    out.pressure_hPa = gas.pressure_hPa(:,ygas);
    if isfield(gas,'lon')
        out.lon = gas.lon(:,ygas);
        out.lat = gas.lat(:,ygas);
    end
    tanstruct = out;
    filename_out = strcat(filein_pre,gasname_out,filein_post);
    savedest = fullfile(matdirectory, filename_out);
    fprintf('\nSaving %s to %s\n', gasname_out, filename_out);
    save(savedest,'tanstruct');
end
%
end


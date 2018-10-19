 function [ smiles_out ] = extract_smiles_data( smiles_in )
%This function reads in a .mat of MLS data and extracts the
%relevant details into one structure.
%
% *INPUT*    
%           mls_in: CELL ARRAY - the MLS data. The .mat file can be created
%           with 'read_mls_data.m'.
%
%
% *OUTPUT*
%           mls_out: STRUCURE - the mls data. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 06/18

%% define some things
%USER DEFINED
% apply_pressure_lim = 0;

smiles_data = smiles_in;
version = strcat('SMILES_', smiles_data{1,1}.band, '_', smiles_data{1,1}.version);
lmls = length(smiles_data);
L2 = [];
prec = [];
lat = [];
lon = [];
time = [];
pressure_hPa = [];
% pressure_hPa = double(smiles_data{1,1}.Pressure);
% lalt = length(pressure_hPa);

%% check to make sure that the files are all for the same gas
gas = cell(1,lmls);
for i = 1:lmls
   gas{i} = smiles_data{1,i}.swath;
end
if length(unique(gas)) > 1
    error('there are measurements for more than one species in the input data')
else
    gas = gas{1};
end

% %% set the flag limits. These numbers are given in the MLS data quality document. Add more cases for another gas that isnt mentioned here
% switch gas
%     case 'ClO'
%         quality_lower_lim = 1.3;
%         convergence_upper_lim = 1.05;
%         p_lower_lim = 147; % in hPa
%         p_upper_lim = 1;
%         % there is also a bias correction that must be done for MLS ClO.
%         % The ASCII file with the correction can be downloaded from the
%         % Aura MLS website
%         load('MLS_ClO_bias_correction.mat'); % loads as 'bias_correction'
%         clo_bias = bias_correction;
%     case 'CO'
%         quality_lower_lim = 1.5;
%         convergence_upper_lim = 1.03;
%         p_lower_lim = 100; % in hPa
%         p_upper_lim = 0.0046;
%     case 'O3'
%         quality_lower_lim = 1.0;
%         convergence_upper_lim = 1.03;
%         p_lower_lim = 100; % in hPa
%         p_upper_lim = 0.02;
%     case 'HOCl'
%         quality_lower_lim = 1.2;
%         convergence_upper_lim = 1.05;
%         p_lower_lim = 10; % in hPa
%         p_upper_lim = 2.2;
% end

%% extract the mls data by measurement day
% if apply_pressure_lim == 1
%     pgood = find(pressure_hPa <= p_lower_lim & pressure_hPa >= p_upper_lim);
% else
%     pgood = find(pressure_hPa);
% end

altitude_km = double(smiles_data{1,1}.Altitude);
pgood = find(altitude_km);
for i = 1:lmls
    if mod(i,100) == 0
        fprintf('past file %i of %i\n',i,lmls)
    end
    valuei = double(smiles_data{1,i}.L2Value(pgood,:));
    preci = double(smiles_data{1,i}.L2Precision(pgood,:));
    lati = double(smiles_data{1,i}.Latitude);
    loni = double(smiles_data{1,i}.Longitude);
    timei = smiles_data{1,i}.Time;
    statusi = smiles_data{1,i}.Status;
    pressurei = double(smiles_data{1,i}.Pressure);
%     qualityi = smiles_data{1,i}.Quality;
%     convergencei = smiles_data{1,i}.Convergence;
    
    if isempty(valuei) == 0
        for j= 1:length(valuei(1,:))
            if statusi(j) == 0 
                L2 = [L2, double(valuei(:,j))];
                prec = [prec, double(preci(:,j))];
                lat = [lat, double(lati(j))];
                lon = [lon, double(loni(j))];
                time = [time,timei(j)];
                pressure_hPa = [pressure_hPa, pressurei(:,j)]; 
            end
        end
    end 
end

%% Change the time to MJD
%SMILES time is seconds from midnight, january 1st, 1958. Obviously.
mjdsmiles = date2mjd(1958,1,1,0,0,time); % this is UTC

%% sort the data by time and put into output structure
[mjdsmiles, I] = sort(mjdsmiles);

smiles_out.source_file = version;
smiles_out.gas = gas;
smiles_out.vmr = L2(:,I);
smiles_out.vmr_error = prec(:,I);
smiles_out.altitude_km = altitude_km;
smiles_out.pressure_hPa = pressure_hPa(:,I);
smiles_out.date_mjd = mjdsmiles;
smiles_out.lat = lat;
smiles_out.lon = lon;

% %% bias correction for ClO
% switch gas
%     case 'ClO'
%        disp('subtracting ClO bias...')
%        %get the indices of the bias correction pressure grid that match the
%        %ordinary grid (current values are 146.78, 100, and 68.13. Indexes = 6, 7, and 8)
%        pressure_hPa_round = round2(pressure_hPa, 0.01); % the bias pressures are rounded to two decimal places
%        ipressure_hPa_bias = find(pressure_hPa_round <= clo_bias.pressure_hPa(1) & pressure_hPa_round >= clo_bias.pressure_hPa(end)); %should be [6,7,8]
%        % interpolate the bias value to the mls measurement latitude and
%        % fill subtract it from the ClO measurements (corresponding to the
%        % above indices).
% %        bias_value_int = nan(length(ipressure_hPa_bias), length(lat));
% %        for n = 1:length(lat)
% %             bias_value_int(:,n) = interp1(clo_bias.lat, clo_bias.value', lat(n));
% %        end
%        bias_value_int = interp1(clo_bias.lat, clo_bias.value', lat');
% %        whos
%        bias_value_int(isnan(bias_value_int)) = 0;
%        mls_out.vmr(ipressure_hPa_bias, :) = mls_out.vmr(ipressure_hPa_bias, :) - bias_value_int';
% end

disp('all done')
%
end


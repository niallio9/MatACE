function [ mls_out ] = extract_mls_data( mls_in )
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

%% define some things
%USER DEFINED
apply_pressure_lim = 1;

mls_data = mls_in;
version = strcat('MLS_', mls_data{1,1}.version);
lmls = length(mls_data);
L2 = [];
prec = [];
lat = [];
lon = [];
time = [];
pressure_hPa = double(mls_data{1,1}.Pressure);
% lalt = length(pressure_hPa);

%% check to make sure that the files are all for the same gas
gas = cell(1,lmls);
for i = 1:lmls
   gas{i} = mls_data{1,i}.swath;
end
if length(unique(gas)) > 1
    error('there are measurements for more than one species in the input data')
else
    gas = gas{1};
end

%% set the flag limits. These numbers are given in the MLS data quality document. Add more cases for another gas that isnt mentioned here
switch gas
    case 'ClO'
        quality_lower_lim = 1.3;
        convergence_upper_lim = 1.05;
        p_lower_lim = 147; % in hPa
        p_upper_lim = 1;
        % there is also a bias correction that must be done for MLS ClO.
        % The ASCII file with the correction can be downloaded from the
        % Aura MLS website
        load('MLS_ClO_bias_correction.mat'); % loads as 'bias_correction'
        clo_bias = bias_correction;
    case 'CO'
        quality_lower_lim = 1.5;
        convergence_upper_lim = 1.03;
        p_lower_lim = 100; % in hPa
        p_upper_lim = 0.0046;
    case 'O3'
        quality_lower_lim = 1.0;
        convergence_upper_lim = 1.03;
        p_lower_lim = 100; % in hPa
        p_upper_lim = 0.02;
    case 'HOCl'
        quality_lower_lim = 1.2;
        convergence_upper_lim = 1.05;
        p_lower_lim = 10; % in hPa
        p_upper_lim = 2.2;
end

%% extract the mls data by measurement day
if apply_pressure_lim == 1
    pgood = find(pressure_hPa <= p_lower_lim & pressure_hPa >= p_upper_lim);
else
    pgood = find(pressure_hPa);
end

pressure_hPa = pressure_hPa(pgood);
for i = 1:lmls
    if mod(i,100) == 0
        fprintf('past file %i of %i\n',i,lmls)
    end
    valuei = mls_data{1,i}.L2gpValue(pgood,:);
    preci = mls_data{1,i}.L2gpPrecision(pgood,:);
    lati = double(mls_data{1,i}.Latitude);
    loni = double(mls_data{1,i}.Longitude);
    timei = mls_data{1,i}.Time;
    statusi = mls_data{1,i}.Status;
    qualityi = mls_data{1,i}.Quality;
    convergencei = mls_data{1,i}.Convergence;
    
    if isempty(valuei) == 0
        for j= 1:length(valuei(1,:))
            if statusi(j) == 0 && qualityi(j) > quality_lower_lim && convergencei(j) < convergence_upper_lim 
                L2 = [L2, double(valuei(:,j))];
                prec = [prec, double(preci(:,j))];
                lat = [lat, double(lati(j))];
                lon = [lon, double(loni(j))];
                time = [time,timei(j)];     
            end
        end
    end 
end

%% Change the time to MJD
%MLS time is seconds from midnight, january 1st, 1993. Obviously.
mjdmls = date2mjd(1993,1,1,0,0,time); % this is UTC

%% sort the data by time and put into output structure
[mjdmls, I] = sort(mjdmls);

mls_out.source_file = version;
mls_out.gas = gas;
mls_out.vmr = L2(:,I);
mls_out.vmr_error = prec(:,I);
mls_out.pressure_hPa = pressure_hPa;
mls_out.date_mjd = mjdmls;
mls_out.lat = lat;
mls_out.lon = lon;

%% bias correction for ClO
switch gas
    case 'ClO'
       disp('subtracting ClO bias...')
       %get the indices of the bias correction pressure grid that match the
       %ordinary grid (current values are 146.78, 100, and 68.13. Indexes = 6, 7, and 8)
       pressure_hPa_round = round2(pressure_hPa, 0.01); % the bias pressures are rounded to two decimal places
       ipressure_hPa_bias = find(pressure_hPa_round <= clo_bias.pressure_hPa(1) & pressure_hPa_round >= clo_bias.pressure_hPa(end)); %should be [6,7,8]
       % interpolate the bias value to the mls measurement latitude and
       % fill subtract it from the ClO measurements (corresponding to the
       % above indices).
%        bias_value_int = nan(length(ipressure_hPa_bias), length(lat));
%        for n = 1:length(lat)
%             bias_value_int(:,n) = interp1(clo_bias.lat, clo_bias.value', lat(n));
%        end
       bias_value_int = interp1(clo_bias.lat, clo_bias.value', lat');
%        whos
       bias_value_int(isnan(bias_value_int)) = 0;
       mls_out.vmr(ipressure_hPa_bias, :) = mls_out.vmr(ipressure_hPa_bias, :) - bias_value_int';
end

disp('all done')
%
end


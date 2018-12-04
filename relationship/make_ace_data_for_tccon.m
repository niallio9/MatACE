function [ ace_tccon_out ] = make_ace_data_for_tccon(  )
%A function to calculate the relationship between two ACE gas measurements.
%The correlation and linear relationship of coincident data points is
%calculated.

% *INPUT*
%           tanstruct1_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           tanstruct1_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with the coincident VMR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 09/2018

%% Define some things
matdirectory = 'C:\Users\ryann\ACE\matdata\';
filein_pre = 'ACE_v3p6_';
filein_post = '.mat';

%% read the data
disp('reading gas data...')
data = load(strcat(matdirectory,filein_pre,'N2O',filein_post));
n2o = apply_ace_flags(data.tanstruct);
data = load(strcat(matdirectory,filein_pre,'CH4',filein_post));
ch4 = apply_ace_flags(data.tanstruct);
data = load(strcat(matdirectory,filein_pre,'HF',filein_post));
hf = apply_ace_flags(data.tanstruct);
data = load(strcat(matdirectory,filein_pre,'T',filein_post));
temp = apply_ace_flags(data.tanstruct);
clear data

%% remove occultations with dodgy pressure values. Just for n2o will work
n2o = filter_ace_pressure(n2o);
% match the occultations of the n2o data and the dmp data
disp('reading DMP data...')
data = load(strcat(matdirectory,filein_pre,'DMP',filein_post));
dmp = remove999_ace_dmp(data.dmpstruct);
clear data
% size(n2o.occultation)
[n2o, dmp] = match_ace_data_dmp(n2o, dmp);
% size(n2o.occultation)

%% match the occultations in each data set
disp('matching the data sets...')
[n2o, ch4] = match_ace_data(n2o, ch4);
[n2o, hf]  = match_ace_data(n2o, hf);
[n2o, temp] = match_ace_data(n2o, temp);
[n2o, ch4] = match_ace_data(n2o, ch4); % gotta do again to make sure any changes are included in ch4
[n2o, hf]  = match_ace_data(n2o, hf); % gotta do again to make sure any changes are included in hf
% size(n2o.occultation)
[n2o, dmp] = match_ace_data_dmp(n2o, dmp); % gotta do again to make sure any changes are included in dmp
% size(n2o.occultation)

% % reduce the altitude range of the datasets to match the dmp data
% disp('subsetting the data to the DMP altitude range...')
% n2o = subset_ace_by_alt(n2o, dmp.altitude_km(1), dmp.altitude_km(end));
% ch4 = subset_ace_by_alt(ch4, dmp.altitude_km(1), dmp.altitude_km(end));
% hf = subset_ace_by_alt(hf, dmp.altitude_km(1), dmp.altitude_km(end));
% temp = subset_ace_by_alt(temp, dmp.altitude_km(1), dmp.altitude_km(end));

%% fix dodgy lat and lon
disp('adding and filtering the location data...')
data = load(strcat(matdirectory,filein_pre,'GLC',filein_post));
n2o = merge_ace_glc(n2o, data.glcstruct);
clear data
% n2o.lat = dmp.lat;
% n2o.lon = dmp.lon;
n2o = filter_ace_bad_lat(n2o);
n2o = filter_ace_bad_lon(n2o);

%% make sure the altitude grids of the dmp and ace data are the same
if ~isequal(n2o.altitude_km(1:length(dmp.altitude_km)), dmp.altitude_km)
    error('the altitude grids of each gas are not equal. Stopping...')
else
    izend = length(dmp.altitude_km);
end

%% start the output structure
disp('making the output structure...')
out.source_file = filein_pre;
out.occultation = n2o.occultation;
out.sr1ss0 = n2o.sr1ss0;
out.date_mjd = n2o.date_mjd;
out.altitude_km = n2o.altitude_km(1:izend);
out.pressure_hPa = n2o.pressure_hPa(1:izend,:);
out.lat = n2o.lat(1:izend,:);
out.lon = n2o.lon(1:izend,:);
out.eql = dmp.eql;
out.spv = dmp.spv;
out.n2o_vmr = n2o.vmr(1:izend,:);
out.n2o_vmr_error = n2o.vmr_error(1:izend,:);
out.ch4_vmr = ch4.vmr(1:izend,:);
out.ch4_vmr_error = ch4.vmr_error(1:izend,:);
out.hf_vmr = hf.vmr(1:izend,:);
out.hf_vmr_error = hf.vmr_error(1:izend,:);
out.temperature_K = temp.vmr(1:izend,:);
% out.temperature_K_error = temp.vmr_error(1:izend,:);


ace_tccon_out = out;
%
disp('Done :)')
%
end


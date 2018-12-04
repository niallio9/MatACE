function [ ] = write_ace_to_netcdf_TCCON( structin, filename_out )
 %A function to create an ACE climatology netCDF file from the information
 %contained in the ACE climstruct .mat file.
%
% *INPUT*
%
%           structin: STRING - the structure containing the gas and DMP
%           information. It can be created using
%           'make_ace_data_for_tccon.m'.
%
%           filename_out: STRING - the name of the .nc file that will be
%           created to hold the data. This can be a directory
%
% *OUTPUT*
%           .nc files of the ACE and DMP information will be written to
%           'filename_out'.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 11/18

%% define some things
data = structin;
filenameout = filename_out;
sdate1950 = datenum(1950,1,1);
timeout = mjd2datenum(data.date_mjd) - sdate1950;

%% change the nans to -999 for bad data
%lat and lon data should contain no nans
data.n2o_vmr(isnan(data.n2o_vmr)) = -999;
data.n2o_vmr_error(isnan(data.n2o_vmr_error)) = -999;
data.ch4_vmr(isnan(data.ch4_vmr)) = -999;
data.ch4_vmr_error(isnan(data.ch4_vmr_error)) = -999;
data.hf_vmr(isnan(data.hf_vmr)) = -999;
data.hf_vmr_error(isnan(data.hf_vmr_error)) = -999;
data.temperature_K(isnan(data.temperature_K)) = -999;
% data.temperature_K_error(isnan(data.temperature_K_error)) = -999;
data.spv(isnan(data.spv)) = -999;
data.eql(isnan(data.eql)) = -999;

%% create and add the data to the file
make_ace_empty_netcdf_TCCON(data,filenameout); % make the empty .nc file (no data)
disp('writing the data to the .nc file')
ncwrite(filenameout, 'occultation', data.occultation);
ncwrite(filenameout, 'sr1ss0', data.sr1ss0);
ncwrite(filenameout, 'time', timeout);
ncwrite(filenameout, 'zlev', data.altitude_km);
ncwrite(filenameout, 'plev', data.pressure_hPa);
ncwrite(filenameout, 'lat', data.lat);
ncwrite(filenameout, 'lon', data.lon);
ncwrite(filenameout, 'eql', data.eql);
ncwrite(filenameout, 'spv', data.spv);
ncwrite(filenameout, 'n2o', data.n2o_vmr);
ncwrite(filenameout, 'n2o_error', data.n2o_vmr_error);
ncwrite(filenameout, 'ch4', data.ch4_vmr);
ncwrite(filenameout, 'ch4_error', data.ch4_vmr_error);
ncwrite(filenameout, 'hf', data.hf_vmr);
ncwrite(filenameout, 'hf_error', data.hf_vmr_error);
ncwrite(filenameout, 'temperature', data.temperature_K);
% ncwrite(filenameout, 'temperature_error', data.temperature_K_error);
%
disp('All done :)')
%
end


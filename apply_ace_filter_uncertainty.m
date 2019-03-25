function [ tanstruct_out_filtered ] = apply_ace_filter_uncertainty(tanstruct_in)
%A function to apply the a filtering whereby only data that lies within 3
%median absolute deviations of the median is kept.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           tanstruct_out_filtered: STRUCTURE - output has the same
%           fields as the input, but filtered according to three median
%           absolute deviations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 05/2018

%% Sort some things
datain = tanstruct_in;
% zlev = length(datain.altitude_km(:,1));
% zlev = length(datain.pressure_hPa(:,1));
out = tanstruct_in; %this structure will be edited later on

%% filter out any data for which the vmr error is larger than 100% of the vmr or less than 0.01%, also where the error is < 0
error_percent = 100 * abs(out.vmr_error)./abs(out.vmr); % put the error as a percentage 
goodcol = find(error_percent < 100 & error_percent > 0.01 & out.vmr_error > 0 );
out = reduce_tanstruct_data_by_index(out, goodcol);
% and where the vmr is > 20ppm, for ozone
if strcmp(out.gas, 'O3')
    goodcol = find(out.vmr < 20e-6 );
    out = reduce_tanstruct_data_by_index(out, goodcol);
end
% whos

%% remove any columns that are only nans now
goodcol = find(nansum(out.vmr) ~= 0);
out = reduce_tanstruct_by_rowindex(out, goodcol);

%%
tanstruct_out_filtered = out;

end


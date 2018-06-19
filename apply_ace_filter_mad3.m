function [ tanstruct_out_filtered ] = apply_ace_filter_mad3(tanstruct_in)
%A function to apply the a filtering whereby only data that lies within 3
%median absolute deviations of the median is kept.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           tanstruct_out_flagged: STRUCTURE - output has the same
%           fields as the input, but filtered according to three median
%           absolute deviations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 05/2018

%% Remove flagged values
datain = tanstruct_in;
% zlev = length(datain.altitude_km(:,1));
zlev = length(datain.pressure_hPa(:,1));
out = tanstruct_in; %this structure will be edited later on

medianvmr = nanmedian(datain.vmr,2); % get the median of the data at each altitude
madvmr = mad(datain.vmr,1,2); % get the median absolute deviation at each altitude
lim_upper = medianvmr + 3.*madvmr;
lim_lower = medianvmr - 3.*madvmr;

%% go through each altitude and remove values outside 3 mads
for n = 1:zlev
    ibad = find(datain.vmr(n,:) > lim_upper(n) | datain.vmr(n,:) < lim_lower(n));
    out.vmr(n,ibad) = nan;
    out.vmr_error(n,ibad) = nan; 
end

%% remove any columns that are only nans now
goodcol = find(nansum(out.vmr) ~= 0);
out = reduce_tanstruct_by_rowindex(out, goodcol);

%%
tanstruct_out_filtered = out;

end


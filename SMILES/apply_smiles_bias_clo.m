function [ smiles_out, bias_out ] = apply_smiles_bias_clo( smiles_in )
%A function to read the SMILES data in .he5 format and output a cell array of
%structure with the information. Use the 'mls_list' input to subset the files
%to a specific year or gas, etc.

% *INPUT*    
%           smiles_list: STRING - the name of (or path to) a text file that
%           contains the names of all of the SMILES files that you would like
%           to read in. The SMILES data generally comes in files that
%           corespond to the day of the year.
%
%           swath: STRING - the name of the gas to which the data
%           corresponds
%
% *OUTPUT*
%           smiles_out: CELL ARRAY - each cell corresponds to one .he5 file. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 05/19

clo = smiles_in; % same fields as for the ACE format
lat_bounds = -90: 10: 90;
clo_out = [];
zunder35 = find(clo.altitude_km < 35);
bias_out = nan(length(zunder35), length(lat_bounds - 1), 12);

for i = 1 : length(lat_bounds) - 1
    fprintf('lat %0.2f - %0.2f\n', lat_bounds(i), lat_bounds(i + 1))
    clo_lati = subset_ace_by_lat_tangent(clo, lat_bounds(i), lat_bounds(i + 1));
    if ~isempty(clo_lati.vmr)
        for j = 1:12
            fprintf('month %i\n', j)
            clo_lati_monthj = subset_ace_by_month(clo_lati, j);
            if ~isempty(clo_lati_monthj.vmr)
                clo_night = subset_ace_by_lst_tangent(clo_lati_monthj, 18, 6); % data between 18:00 and 6:00
                bias_lati_monthj = nanmean(clo_night.vmr(zunder35, :), 2);
                bias_lati_monthj(isnan(bias_lati_monthj)) = 0;
                clo_lati_monthj.vmr(zunder35, :) = clo_lati_monthj.vmr(zunder35, :) - bias_lati_monthj; % subtract the nighttime mean
                clo_out = merge_ace_data(clo_out, clo_lati_monthj);
                bias_out(:, i, j) = bias_lati_monthj;
            end
        end
    end 
end
smiles_out = clo_out;
disp('done')
%
end


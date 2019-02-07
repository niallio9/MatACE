function [ mlsstruct_out ] = merge_mls_dmp( mlsstruct_in, dmpstruct_in, save_appendix )
%A function to match the occultations of the ACE data and dmp structures.
%This is done using the intersection of occultation numbers.

% *INPUT*
%           mlsstruct_in: STRUCTURE - a .MAT structure containing MLS
%           data and metadata. It is usually created using
%           'read_and_extract_mls_data'.
%
%           dmpstruct_in: STRUCTURE - a .MAT structure containing MLS DMP
%           data and metadata. It is usually created using
%           'read_and_extract_mls_data'.
%
% *OUTPUT*
%           mlsstruct_out: STRUCTURE - output has the same
%           fields as the input, but with the data that coincides with the
%           dmp info.
%
%           dmpstruct_out: STRUCTURE - output has the same
%           fields as the input, but with the data that coincides with the
%           gas info.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 12/2018

%% Define some things
if nargin < 3
    tosave = 0;
    save_appendix = '';
else
    tosave = 1;
    if ~isempty(save_appendix)
        save_appendix = strcat('_',save_appendix);
    end
end
interptype = 'pchip';
gas = mlsstruct_in;
dmp = dmpstruct_in;
clear mlsstruct_in
gas.spv = nan(size(gas.vmr));

[ ~, ~, ygas, ydmp ] = match_mls_data_dmp( gas, dmp ); % these are the repective indices of the data that match
disp('interpolating the sPV values to the pressure levels of the gas data...')
gas.spv(:,ygas) = interp1(dmp.pressure_hPa, dmp.spv(:,ydmp), gas.pressure_hPa, interptype, nan); % interpolate the spv to the gas pressure grid
disp('done')

mlsstruct_out = gas;

%% save the data
if tosave == 1
    mlsstruct = mlsstruct_out; % for the naming of the output variable
    savedest = fullfile(pwd, strcat('MLS_v4p2_', mlsstruct_out.gas,'_dmp', save_appendix, '.mat'));
    fprintf('saving data to %s\n', savedest);
    save(savedest,'mlsstruct','-v7.3');
    fprintf('done\n')
end

disp('all done :)')
%
end


function [ mlsstruct_out ] = merge_mls_dmp_by_year( mlsstruct_in, years_in, to_save )
%A function to match the occultations of the ACE data and dmp structures.
%This is done using the intersection of occultation numbers.

% *INPUT*
%           mlsstruct_in: STRUCTURE - a .MAT structure containing MLS
%           data and metadata. It is usually created using
%           'read_and_extract_mls_data'.
%
%           years_in: STRUCTURE - the years of the data that you want to
%           find matches for.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
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
else
    tosave = to_save;
end
% dmp_dir = 'C:\Users\ryann\MLS\dmpdata\';
dmp_dir = '/net/deluge/pb_1/users/nryan/MLS/dmpdata/';
interptype = 'pchip';
if ~isdir(dmp_dir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n', dmp_dir)
    error('The DMP directory couldn''t be found')
end
gas = mlsstruct_in;
clear mlsstruct_in
gas.spv = nan(size(gas.vmr));
yearsin = years_in;
lyears = length(yearsin);

%% get the matching spv data
for i = 1:lyears
    fprintf('\nFor %i...\n',yearsin(i))
    load(strcat(dmp_dir, 'MLS_v4p2_DMP_', num2str(yearsin(i)), '.mat'))
    dmp = mlsstruct; % 'mlsstruct' will be the name of the loaded structure
    [ ~, ~, ygas, ydmp ] = match_mls_data_dmp( gas, dmp ); % these are the repective indices of the data that match
    disp('interpolating the sPV values to the pressure levels of the gas data...')
    gas.spv(:,ygas) = interp1(dmp.pressure_hPa, dmp.spv(:,ydmp), gas.pressure_hPa, interptype, nan); % interpolate the spv to the gas pressure grid
    disp('done')
end
clear mlsstruct dmp
mlsstruct_out = gas;

%% save the data
if tosave == 1
    mlsstruct = mlsstruct_out; % for the naming of the output variable
    savedest = fullfile(pwd, strcat('MLS_v4p2_', mlsstruct_out.gas,'_dmp',num2str(yearsin(1)), num2str(yearsin(end)),'.mat'));
    fprintf('saving data to %s\n', savedest);
    save(savedest,'mlsstruct','-v7.3');
    fprintf('done\n')
end

disp('all done :)')
%
end


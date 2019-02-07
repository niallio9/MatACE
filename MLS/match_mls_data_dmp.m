function [ mlsstruct_out, dmpstruct_out, ygas, ydmp ] = match_mls_data_dmp( mlsstruct_in, dmpstruct_in )
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
% dmp_dir = 'C:\Users\ryann\MLS\dmpdata\';
% if ~isdir(dmp_dir)
%     fprintf('\nIt doesn''t look like ''%s'' exists...\n', dmp_dir)
%     error('The DMP directory couldn''t be found')
% end
gas = mlsstruct_in;
dmp = dmpstruct_in;

vdate_gas = nan(length(gas.date_mjd), 8);
vdate_gas(:,1:6) = datevec(mjd2datenum(gas.date_mjd));
vdate_gas(:,7) = round2(gas.lat,1);
vdate_gas(:,8) = round2(gas.lon,1);
vdate_gas(:,6) = [];

vdate_dmp = nan(length(dmp.date_mjd), 8);
vdate_dmp(:,1:6) = datevec(mjd2datenum(dmp.date_mjd));
vdate_dmp(:,7) = round2(dmp.lat,1);
vdate_dmp(:,8) = round2(dmp.lon,1);
vdate_dmp(:,6) = [];

% check to see if they already match
if isequal(vdate_gas(:,1:4), vdate_dmp(:,1:4))
    disp('the occultation numbers already match')
    mlsstruct_out = gas;
    dmpstruct_out = dmp;
else
    fprintf('\nsubsetting the mls data and DMPs by matching times down to the minute, and the lat and lon\n')
    
    % Check which ones match
    [~,ygas,ydmp] = intersect(vdate_gas, vdate_dmp, 'rows'); % the indices of where the year/month/date/hour/minute match
    
    %% reduce the sizes of the variables to only include the ones with coincident observations
    mlsstruct_out = gas;
    clear mlsstruct
    mlsstruct_out.vmr = mlsstruct_out.vmr(:,ygas);
    mlsstruct_out.vmr_error = mlsstruct_out.vmr_error(:,ygas);
    mlsstruct_out.date_mjd = mlsstruct_out.date_mjd(ygas);
    mlsstruct_out.lat = mlsstruct_out.lat(ygas);
    mlsstruct_out.lon = mlsstruct_out.lon(ygas);
    
    dmpstruct_out = dmp;
    clear dmp
    dmpstruct_out.spv = dmpstruct_out.spv(:,ydmp);
    dmpstruct_out.date_mjd = dmpstruct_out.date_mjd(ydmp);
    dmpstruct_out.lat = dmpstruct_out.lat(ydmp);
    dmpstruct_out.lon = dmpstruct_out.lon(ydmp);
    fprintf('Done\n')
end
%
end


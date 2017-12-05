function [ tanstruct_out, dmpstruct_out ] = match_ace_data_dmp( tanstruct_in, dmpstruct_in )
%A function to match the occultations of the ACE data and dmp structures.
%This is done using the intersection of occultation numbers.

% *INPUT*
%           tanstruct_in: STRUCTURE - a .MAT structure containing ACE
%           data and metadata. It is usually created using
%           'read_ace_ncdata', 'read_ace_ncdata_for_mat'.
%
%           dmpstruct_in: STRUCTURE - containins the ACE dmp data. It is
%           usually created using 'read_ace_dmp' or 'read_ace_dmp_for_mat'.
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
% NJR - 10/2017

%% Define some things
gas = tanstruct_in;
dmp = dmpstruct_in;
lgas = length(gas.occultation);
ldmp = length(dmp.occultation);

%% check to see if they already match
if isequal(gas.occultation, dmp.occultation)
    %disp('the occultation numbers already match')
    tanstruct_out = gas;
    dmpstruct_out = dmp;
else
    fprintf('\nsubsetting the ace data and DMPs by matching occultation numbers\n')
    %% Get the unique codes that identify the occultations/orbit.
    gasorbit(1,1:lgas) = gas.occultation;
    gasorbit(2,1:lgas) = gas.sr1ss0;
    dmporbit(1,1:ldmp) = dmp.occultation;
    dmporbit(2,1:ldmp) = dmp.sr1ss0;
    
    %% Check which ones match
    [~,ygas,ydmp] = intersect(gasorbit',dmporbit','rows'); % the indices of where the orbits/occultations match
    
    %% reduce the sizes of the variables to only include the ones with coincident orbit/occulations
    gasout = reduce_tanstruct_by_rowindex(gas,ygas);
    dmpout = reduce_dmpstruct_by_rowindex(dmp,ydmp);
    
    %% out
    tanstruct_out = gasout;
    dmpstruct_out = dmpout;
    fprintf('Done\n')
end
%
end


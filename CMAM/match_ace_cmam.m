function [ tanstruct_out, cmam_out ] = match_ace_cmam( tanstruct_in, cmam_in )
%A function to match the sizes of ACE data-type structures. This is
%done using the intersection of the occultation numbers.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           cmam_in: STRUCTURE - contains the gas specific CMAM data that
%           has aready been sampled to the location of the ACE data. The
%           structure created by 'sample_cmam_for_ace.m' can also be used
%           to create the structure.
%
% *OUTPUT*
%           tanstruct1_out: STRUCTURE - output has the same fields as the
%           input, but with the data that coincides with the other input.

%           tanstruct2_out: STRUCTURE - same as above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 10/2017

%% Define some things
gas1 = tanstruct_in;
gas2 = cmam_in;
lgas1 = length(gas1.occultation);
lgas2 = length(gas2.occultation);

%% check to see if they already match
if isequal(gas1.occultation, gas2.occultation)
    %disp('the occultation numbers already match')
    tanstruct_out = gas1;
    cmam_out = gas2;
else
    fprintf('\nSubsetting the ace data by matching occultation numbers...\n')
    %% Get the unique codes that identify the occultations/orbit.
    gas1orbit(1,1:lgas1) = gas1.occultation;
    gas1orbit(2,1:lgas1) = gas1.sr1ss0;
    gas2orbit(1,1:lgas2) = gas2.occultation;
    gas2orbit(2,1:lgas2) = gas2.sr1ss0;
    
    %% Check which ones match
    [~,ygas1,ygas2] = intersect(gas1orbit',gas2orbit','rows'); % the indices of where the orbits/occultations match
    
    %% reduce the sizes of the variables to only include the ones with coincident orbit/occulations
    gas1out = reduce_tanstruct_by_rowindex(gas1,ygas1);
    gas2out = reduce_cmam_by_rowindex(gas2,ygas2);
    
    %% out
    tanstruct_out = gas1out;
    cmam_out = gas2out;
    fprintf('Done\n')
end
%
end


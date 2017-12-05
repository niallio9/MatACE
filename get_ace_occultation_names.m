function [ occultation_names ] = get_ace_occultation_names( tanstruct_in )
%A function to make the occultation names for a set of ACE measurements.
%The naming convention is used for a few other things like the geolocation
%files.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. Structures containing information
%           about the ACE DMPs (dmpstruct) and ACE GLCs can also be input
%           here.
%
%
% *OUTPUT*
%           occultation_names: STRING - contains the occultation
%           names that correspond to the data in the input  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 11/17

%% Define some things
gas = tanstruct_in;
lorbit = length(gas.occultation);
namesout = strings(1,lorbit);

%% look for the right file according to the info in tanstruct
orbits = gas.occultation;
srss = gas.sr1ss0;
for i = 1:lorbit
    orbits(i);
    if srss(i) == 1
        namesout{i} = strcat('sr',num2str(orbits(i)));
    elseif srss(i) == 3
        namesout{i} = strcat('sr',num2str(orbits(i)),'a');
    elseif srss(i) == 0
        namesout{i} = strcat('ss',num2str(orbits(i)));
    elseif srss(i) == 2
        namesout{i} = strcat('ss',num2str(orbits(i)),'a');
    else
        sprintf('there''s no info for ')
        namesout{i} = nan;
    end
end
occultation_names = namesout;
%
end


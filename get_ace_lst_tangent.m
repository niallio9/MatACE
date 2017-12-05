function [ lst_out ] = get_ace_lst_tangent( tanstruct_in )
%A function to calculate the local solar time (LST) of ace measurements
%using the ace measurement times and the tangent longitudes.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           lst_out: VECTOR - a vector of the same size as the
%           tanstruct_in.date_mjd variable, containing the local solar
%           times of each ace measurement.

%% Define some things
gas = tanstruct_in;

%% create the lst vector

lst = mjd2lst(gas.date_mjd, gas.lon_tangent);

lst_out = lst;
%
end


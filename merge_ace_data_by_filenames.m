function [ tanstruct_out ] = merge_ace_data_by_filenames( filenames_in )
%A function to merge two ace files so that the data is a combination of the
%data in both files. Duplicate occultations are removed

% *INPUT*
%           filenames_in: ARRAY OF CELLS - contains the gas specific ACE data.
%           This structure 
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output is a structure with the data
%           from each of the input files, combined.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 12/2018

%% Define some things

%% Load the first file
fprintf('\nloading %s...\n', filenames_in{1})
load(filenames_in{1}) % variable is called mlsstruct_acesample
disp('done')
tanstruct_out = mlsstruct_acesample;

for i = 2:length(filenames_in)
    fprintf('\nloading %s...\n', filenames_in{i})
    load(filenames_in{i}) % variable is called mlsstruct_acesample
    disp('done')
    
    tanstruct_out = merge_ace_data(tanstruct_out, mlsstruct_acesample);  
end

%%
disp('all done')
%
end


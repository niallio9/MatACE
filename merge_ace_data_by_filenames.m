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
struct_in = load(filenames_in{1}); % variable is called mlsstruct_acesample
struct_name_1 = fieldnames(struct_in);
disp('done')
tanstruct_out = struct_in.(char(struct_name_1)); % set up the first structure to be merged

for i = 2:length(filenames_in)
    fprintf('\nloading %s...\n', filenames_in{i})
    struct_in = load(filenames_in{i});
    struct_name_i = fieldnames(struct_in);
    if ~strcmp(char(struct_name_i), char(struct_name_1))
        warning('the input files don''t have the same structure names')
        fprintf('file 1: %s. file_%i: %s', char(struct_name_1), i, char(struct_name_i))
    end
    tanstruct_i = struct_in.(char(struct_name_i));
    disp('done')
    tanstruct_out = merge_ace_data(tanstruct_out, tanstruct_i);
end

%%
disp('all done')
%
end


function [ ] = read_and_extract_mls_data_by_year( years_in, swath )
%A function to read and extract the MLS data and save it to files by year,
%to the directory specified below.
%It is assumed that you are in the directory of the .he5 files. 

% *INPUT*
%           years_in: STRUCTURE - the years of the data that you want to
%           read.
%
%           swath: STRING - the name of the gas to which the data
%           corresponds
%
% *OUTPUT*
%           The save directory is defined in 'read_and_extract_mls_data.m',
%           which is called below.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 03/2019

%% Define some things
% out_dir = 'C:\Users\ryann\MLS\matdata\';
% if ~isdir(out_dir)
%     fprintf('\nIt doesn''t look like ''%s'' exists...\n', dmp_dir)
%     error('The DMP directory couldn''t be found')
% end

yearsin = years_in;
lyears = length(yearsin);
savename_pre = 'MLS_v4p2_';
savename_post = '.mat';

for i = 1:lyears
    fprintf('\nDoing %i...', yearsin(i))
    fprintf('\nmaking a list for the %i files...\n', yearsin(i))
    listname = sprintf('list_%i.txt', years_in(i));
    make_list(sprintf('*%i*.he5', yearsin(i)), listname);
    
    if strcmp(swath, 'ScaledPV')
        savename = strcat(savename_pre, 'DMP', sprintf('_%i', yearsin(i)), savename_post);
    else
        savename = strcat(savename_pre, swath, sprintf('_%i', yearsin(i)), savename_post);
    end
    read_and_extract_mls_data( listname, swath, savename);
end
%
end


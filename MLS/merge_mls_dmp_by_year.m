function [ ] = merge_mls_dmp_by_year( gas_name, years_in, save_appendix )
%A function to match the occultations of the ACE data and dmp structures.
%This is done using the intersection of occultation numbers.

% *INPUT*
%           mgas_name: STRUCTURE - the name of the gas about which the MLS
%           file holds the data.
%
%           years_in: STRUCTURE - the years of the data that you want to
%           find matches for.
%
%           save_appendix: STRING - an appendix on the name of the saved
%           output file. Can be left blank.
%
% *OUTPUT*
%           Merged files are saved to the current directory with
%           'merge_mls_dmp.m'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 12/2018

%% Define some things
dmp_dir = 'C:\Users\ryann\MLS\dmpdata\';
gas_dir = 'C:\Users\ryann\MLS\matdata\';
if ~isdir(dmp_dir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n', dmp_dir)
    error('The DMP directory couldn''t be found')
end
if ~isdir(gas_dir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n', gas_dir)
    error('The gas directory couldn''t be found')
end
yearsin = years_in;
lyears = length(yearsin);
filename_pre = 'MLS_v4p2_';
filename_post = '.mat';
if nargin < 3
    save_appendix = '';
else
    if ~isempty(save_appendix)
        save_appendix = strcat('_',save_appendix);
    end
end

for i = 1:lyears
    fprintf('\nDoing %i...\n', yearsin(i))
    saveappendix = strcat(num2str(yearsin(i)), save_appendix ); % put the year in the save name
    gasfile = fullfile(gas_dir, strcat(filename_pre, gas_name, '_', num2str(yearsin(i)), filename_post));
    if exist(gasfile,'file') ~= 2
        error('cannot find %s', gasfile)
    end
    dmpfile = fullfile(dmp_dir, strcat(filename_pre, 'DMP', '_', num2str(yearsin(i)), filename_post));
    if exist(dmpfile,'file') ~= 2
        error('cannot find %s', dmpfile)
    end
    load(dmpfile);
    dmpstruct = mlsstruct; % the variable loaded will be mlsstruct
    load(gasfile); % the variable loaded will be mlsstruct
    merge_mls_dmp( mlsstruct, dmpstruct, saveappendix );
end
%
end


function [ ] = sample_and_scale_mls_for_ace_by_year(tanstruct_in, years_in, save_appendix)
%A function to sample MLS data according to the time/lat/lon/alt of ACE
%measurements. Data from a chemical box model is to scale the data.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           save_appendix: STRING - an appendix on the name of the saved
%           output file. Can be left blank.
%
% *OUTPUT*
%           cmam_sample: STRUCTURE - with the data that has been sampled
%           from CMAM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 12/18

%% Define some things
acedir = 'C:\Users\ryann\ACE\matdata\';
mlsdir = 'C:\Users\ryann\MLS\matdata\';
pratfile = 'ACE_v3p6_pratmo_ClO_LST_20042018.mat'; % this is gas-specific !!!
% pratfile = 'ACE_v3p6_pratmo_HOCl_LST_20042018.mat'; % this is gas-specific !!!
mlsfile_pre = 'MLS_v4p2_ClO_dmp_';
% mlsfile_pre = 'MLS_v4p2_HOCl_dmp_';
mlsfile_post = '.mat';
% use_maxmin = 1; % you set this to 1 use the max and min output value in the 'sample_and_scale_mls_for_ace.m' function below. Set to 0 otherwise

%% Check the directories
if ~isdir(acedir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n', acedir)
    error('The ACE directory couldn''t be found')
end
if ~isdir(mlsdir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n', mlsdir)
    error('The MLS directory couldn''t be found')
end
pratfile = fullfile(acedir, 'pratmo', pratfile);
if exist(pratfile,'file') ~= 2
    error('cannot find %s', pratfile)
end

%%
yearsin = years_in;
lyears = length(yearsin);
if nargin < 3
    save_appendix = '';
end
%% main loop
disp('loading the pratmo data...')
load(pratfile); % variable is called pratstruct
disp('done')
for i = 1:lyears
    fprintf('\nDoing %i...\n', yearsin(i))
    if isempty(save_appendix)
        saveappendix = num2str(yearsin(i));
    else
        saveappendix = strcat(save_appendix, '_', num2str(yearsin(i))); % put the year in the save name
    end
    
    %% subset ace data to the year
    disp('subsetting the ace data by year...')
    acestruct_i = subset_ace_by_year(tanstruct_in, yearsin(i));
    disp('done')
    
    %% read the mls data
    mlsfile = fullfile(mlsdir, strcat(mlsfile_pre, num2str(yearsin(i)), mlsfile_post));
    fprintf('\nloading %s...\n', mlsfile)
    if exist(mlsfile,'file') ~= 2
        error('cannot find %s', mlsfile)
    end
    load(fullfile(mlsdir, strcat(mlsfile_pre, num2str(yearsin(i)), mlsfile_post))); % variable is called mlsstruct
    mlsstruct_i = convert_mls_to_ace_format(mlsstruct); % convert here so can use ace-related functions on the structure
    
    %% get the max and min values for the year after filtering out outliers
    disp('getting the max and min of the mls data without outliers')
    mls_minmax.lat_bounds = -90:5:90;
    [ mls_minmax.max_val, mls_minmax.min_val, mls_minmax.lat_bins ] = get_ace_maxmin_by_lat_tangent( mlsstruct_i, mls_minmax.lat_bounds, 'fraction' );
    disp('done')
%     mls_minmax = [];
    
    %% subset both structures by am/pm and run the sampling code twice
    disp('splitting ace and mls data by am/pm...')
    [acestruct_i_am, acestruct_i_pm] = split_ace_by_lst_tangent(acestruct_i);
    clear acestruct_i 
    [mlsstruct_i_am, mlsstruct_i_pm] = split_ace_by_lst_tangent(mlsstruct_i);
    clear mlsstruct_i
    disp('done')
    disp('applying sampling function to the am data...')
    saveappendix_am = strcat('am_',saveappendix);
    sample_and_scale_mls_for_ace(mlsstruct_i_am, acestruct_i_am, pratstruct, saveappendix_am, mls_minmax);
    disp('applying sampling function to the pm data...')
    saveappendix_pm = strcat('pm_',saveappendix);
    sample_and_scale_mls_for_ace(mlsstruct_i_pm, acestruct_i_pm, pratstruct, saveappendix_pm, mls_minmax);
%     mlsstruct_acesample_i = sample_and_scale_mls_for_ace(mlsstruct, acestruct_i, pratstruct, saveappendix);
%     if i == 1
%         mlsstruct_acesample_out = mlsstruct_acesample_i;
%     else
%         mlsstruct_acesample_out = merge_ace_data(mlsstruct_acesample_out, mlsstruct_acesample_i);
%     end
end
%
end

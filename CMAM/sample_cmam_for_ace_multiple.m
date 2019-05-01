function [cmam_acesample] = sample_cmam_for_ace_multiple(tanstruct_in, gasname_in)
%A function to sample cmam data according to the time/lat/lon/alt of ACE
%measurements. The two closest times are sampled, mainly to have a choice
%of start-time when using chemical box model to scale the data at a later
%date.

% *INPUT*    
%
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           gasname_in: STRING - the name of the gas for the file that holds
%           the cmam data. The spelling is all lower case.
%
% *OUTPUT*
%           cmam_sample: STRUCTURE - with the data that has been sampled
%           from CMAM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 07/18
%% Define some things
cmam_dir = 'D:\CMAM\';
ace = tanstruct_in;
gasname_cmam = gasname_in;
% the files are named like:
% 'vmrclo_6hrChem_CMAM_CMAM30-SD_r1i1p1_2004010100-2004123118.nc'
filein_pre = sprintf('vmr%s_6hrChem_CMAM_CMAM30-SD_r1i1p1_',gasname_cmam);
% get the unique years for the data
vdates = datevec(mjd2datenum(ace.date_mjd));
years_unique = unique(vdates(:,1));
nyears = length(years_unique);
cmam_sample_out = [];

for n = 1:nyears % the cmam files are split by year
    % subset the ace data to the year in question
    ace_n = subset_ace_by_year(ace,years_unique(n));
    % get the cmam file for that year
    filename_n = sprintf('%s%i010100-%i123118.nc', filein_pre, years_unique(n), years_unique(n)); % 'vmrclo_6hrChem_CMAM_CMAM30-SD_r1i1p1_2004010100-2004123118.nc'
    filein = strcat(cmam_dir,filename_n); % full path to the file
    temp_dir = dir(filein); % to check whether the file exists
    if(~isempty(temp_dir)) % if the file exists
        fprintf('\nsampling cmam data for %i\n', years_unique(n));
        tic
        cmam_sample_n = sample_cmam_for_ace(filein, ace_n);
        toc
        if n > 1
            disp('merging with previous data...')
        end
        cmam_sample_out = merge_ace_data(cmam_sample_out, cmam_sample_n); 
    else
        fprintf('Error: I can''t find %s\n',filein); % if the file doesn't exist
        fprintf('moving on...\n')
    end
end
%% save the output file
cmam_acesample = cmam_sample_out;
savename = sprintf('cmam_sample_ace_%s_%i%i', gasname_cmam, years_unique(1), years_unique(end));
savedest = strcat(cmam_dir,savename);
fprintf('saving output file %s...', savename)
save(savedest,'cmam_acesample')
disp('done')
%
end

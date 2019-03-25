function [ ] = make_ace_gas_vmrs_with_pratmo_by_year(tanstruct_o3_in, tanstruct_T_in, years_in, save_appendix, gasnames_in)
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

%%
yearsin = years_in;
lyears = length(yearsin);
% if nargin < 4
%     save_appendix = '';
% end
%% main loop
for i = 1:lyears
    fprintf('\nDoing %i...\n', yearsin(i))
    if isempty(save_appendix)
        save_appendix_yeari = num2str(yearsin(i));
    else
        save_appendix_yeari = strcat(save_appendix, '_', num2str(yearsin(i))); % put the year in the save name
    end
    
    %% subset ace data to the year
    disp('subsetting the ace data by year...')
    tanstruct_o3_yeari = subset_ace_by_year(tanstruct_o3_in, yearsin(i));
    tanstruct_T_yeari = subset_ace_by_year(tanstruct_T_in, yearsin(i));
    disp('done')
    
    %% run the box model code
    make_ace_gas_vmrs_with_pratmo( tanstruct_o3_yeari, tanstruct_T_yeari, save_appendix_yeari, gasnames_in )
end
%
end

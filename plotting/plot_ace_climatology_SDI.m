function [ vmrzon_dif, vmrzonvar_dif, obscount_dif, vmrzon_dif_percent, mean_vmr ] = plot_ace_climatology_SDI( filename, gasname_in )
%A funcion to compare the climatology made by Jaho, with the current
%version of the climatology. The assumption is that the two versions are
%made on the same latitude and altitude grid, and ar for the same gas.

% *INPUT*
%           filename_new: string - the .mat file that contains the
%           new version of the climatology data.
%
%           filename_old: string - the netcdf file that contains the
%           old version of the climatology data.
%
% *OUTPUT*
%           vmrzon_dif: array - the differnce of the new and old zonal
%           vmrs.
%
%           vmrzonvar_dif: array - the differnce of the new and old
%           standard deviation of the zonal vmrs.
%
%           obscount_dif: array - the differnce of the new and old
%           observation counts.
%
%           vmrzon_dif_percent: array - the percentage differnce of the new
%           and old zonal vmrs. The new vmrs are the denominator.
%
%           mean_vmr: vector - a mean altitude profile of the new version of
%           the data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *INPUT*
%           climstruct_in: STRUCTURE - contains the gas specific ACE
%           climatology data. This structure can be created with
%           'make_ace_climatology.m' or with
%           'make_ace_climatology_multiple.m'. 
%
% *OUTPUT*
%           oldclim_in: netcdf - contains the gas specific climatology data
%           of the previous version of the climatology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define some things
gasname = gasname_in;
vmrzon = ncread(filename,gasname);
vmrzon = vmrzon(:,6:33,:);
vmrzon(vmrzon == -999) = nan;
% vmrzonvar = ncread(filename_new,strcat(gasname,'_STD')); % is the variance value in the new data
% vmrzonvar = vmrzonvar(:,6:33,:);
% vmrzonvar(vmrzonvar == -999) = nan;
lat = ncread(filename_new,'lat');
p = ncread(filename_new,'plev');
p = p(6:33);
% obscount = ncread(filename_new,strcat(gasname,'_NR'));
% obscount = obscount(:,6:33,:);
% obscount(obscount == -999) = nan;
if strcmp(gasname(end-1:end),'am') || strcmp(gasname(end-1:end),'pm')
    gasname_old = gasname(1:end-3);
else
    gasname_old = gasname;
end

%% compare the fields
vmrzon_dif = vmrzon - vmrzon_old;
vmrzonvar_dif = vmrzonvar - vmrzonvar_old;
obscount_dif = obscount - obscount_old;
vmrzon_dif_percent = 100*vmrzon_dif./vmrzon;
mean_vmr = nanmean(vmrzon,3);

end


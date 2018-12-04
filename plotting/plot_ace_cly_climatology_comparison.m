function [ test, test1, test2 ] = plot_ace_cly_climatology_comparison( do_plot )
%A funcion to compare the climatology made by Jaho, with the current
%version of the climatology. The assumption is that the two versions are
%made on the same latitude and altitude grid, and ar for the same gas.

% *INPUT*
%           gasname: STRING - the name of the gas for which you want to
%           compare the climatologies.
%
%           filename_oldclim: STRING - the netcdf file that contains the
%           old version of the climatology data.
%
% *OUTPUT*
%           vmrzon_dif: CELL OF ARRAYS - the differnce of the new and old
%           zonal vmrs for each month of the year.
%
%           vmrzonvar_dif: CELL OF ARRAYS - the differnce of the new and
%           old standard deviation of the zonal vmrs for each month of the
%           year.
%
%           obscount_dif: CELL OF ARRAYS - the differnce of the new and old
%           observation counts for each month of the year.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define which species you will be plotting

file_pre = 'ACEFTS_CLIM_v3_lat_'; % ACEFTS_CLIM_v3_lat_O3_DJF.mat
file_post = '_20042010.mat';

% dataX = {clo, hocl, hcl, clono2, Cly}; the Cly name is a dummy for later
gasnames1 = {'ClOmlspratlatnegfixampm','HOClmls_sap','HCl','ClONO2'};
gasnames2 = {'ClOcmam','HOClcmam','HClcmam','ClONO2cmam'};

%% define some things
home_linux = '/home/niall/Dropbox/climatology/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\'; %#ok<NASGU>
% newclim_dir = strcat(home_windows,'climdata_v3p5_nr\');
clim_dir = 'C:\Users\ryann\ACE\climdata_testing\time_matched_climatology\';
% newclim_dir = 'C:\Users\ryann\MLS\climdata\';
% newclim_dir = 'C:\Users\ryann\ACE\MAESTRO\climdata\';

if nargin > 1
    yplot = do_plot;
else
    yplot = 1;
end

vmrzon1 = nan(48,36,length(gasnames1) + 1);
vmrzon2 = nan(48,36,length(gasnames1) + 1);

%% loop through the months and fill in the cells of arrays
for i = 1:length(gasnames1)% for each gas
    filei = strcat( clim_dir, gasnames1{i},'/', file_pre, gasnames1{i}, file_post);
    if exist(filei,'file') ~= 2
        fprintf('There is no file for %i_%0.2i. Moving on...\n', yearin, months(i))
        vmrzon1(:,:,i) = nan(48,36); % this is hardcoded here, for now.
    else
        clim = load(filei); clim = clim.climstruct; % the variable is called climstruct in the new data
        vmrzon1(:,:,i) = clim.vmr_zonal;
    end
end
for i = 1:length(gasnames2)% for each gas
    filei = strcat( clim_dir, gasnames2{i},'/', file_pre, gasnames2{i}, file_post);
    if exist(filei,'file') ~= 2
        fprintf('There is no file for %i_%0.2i. Moving on...\n', yearin, months(i))
        vmrzon2(:,:,i) = nan(48,36); % this is hardcoded here, for now.
    else
        clim = load(filei); clim = clim.climstruct; % the variable is called climstruct in the new data
        vmrzon2(:,:,i) = clim.vmr_zonal;
    end
end
pace = clim.pressure_hPa;
zace = clim.altitude_km_mean;
izace_30up = find(zace>30);
latace = clim.lat;

%% replace nans with zeros in the vmr data so that we can ignore points
% get locations of all-nan profiles
vmrzon_hocl_nan2zero = vmrzon1(:,:,2);
vmrzon_clono2_nan2zero = vmrzon1(:,:,4);

Jnanprofile_hocl = nansum(vmrzon_hocl_nan2zero,1) == 0; % index of all-nan columns
Jnanprofile_clono2 = nansum(vmrzon_clono2_nan2zero,1) == 0;
vmrzon_hocl_nan2zero(:,Jnanprofile_hocl) = 999;
vmrzon_clono2_nan2zero(:,Jnanprofile_clono2) = 999;
%hocl: assume missing values values are negligible
vmrzon_hocl_nan2zero(isnan(vmrzon_hocl_nan2zero)) = 0;
%clono2: upper scaled a priori sucks so we cant use it. assume missing
%values above 30km are negligible.
dummy = vmrzon_clono2_nan2zero(izace_30up,:,:);
dummy(isnan(dummy)) = 0.001e-9;
vmrzon_clono2_nan2zero(izace_30up,:,:) = dummy;
% restore the 999 values to nans.
vmrzon_hocl_nan2zero(vmrzon_hocl_nan2zero == 999) = nan;
vmrzon_clono2_nan2zero(vmrzon_clono2_nan2zero == 999) = nan;
% add the altered data back into the vmrzon array
vmrzon1_cly = vmrzon1;
vmrzon1_cly(:,:,2) = vmrzon_hocl_nan2zero;
vmrzon1_cly(:,:,4) = vmrzon_clono2_nan2zero;

%clean up the data
vmrzon1(vmrzon1 < 0.001e-9) = 0.001e-9;
vmrzon2(vmrzon2 < 0.001e-9) = 0.001e-9;
vmrzon1_cly(vmrzon1_cly < 0.001e-9) = 0.001e-9;

% make Cly
vmrzon1(:,:,5) = sum(vmrzon1_cly(:,:,1:4),3);
vmrzon2(:,:,5) = sum(vmrzon2(:,:,1:4),3);

%get the differences in the data
vmrzondif = vmrzon1 - vmrzon2;
meanvmrzon = (vmrzon1+vmrzon2)/2;
vmrzondifp = 100 * vmrzondif ./ meanvmrzon;
max(vmrzondif(:));
max(vmrzondifp(:))
% vmrzon_cly = squeeze(vmrzon1(:,:,1)) + squeeze(vmrzon1(:,:,2)) + squeeze(vmrzon1(:,:,1)) + squeeze(vmrzon1(:,:,1));

%% Make the plots if you want
vmrzon1 = vmrzon1.*1e9;
vmrzon2 = vmrzon2.*1e9;
vmrzondif = vmrzondif.*1e9;

c = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.5 1 2 3 4];
cgrid = log(c);

cminus = [0.001 0.01 0.1 1];
% cminus = logspace(0.02,max(vmrzondif(:)),6);
% cminus = cminus(1:end-3);
cminus = [fliplr(cminus * -1),cminus];
cminusgrid = log(abs(cminus));
mincminusgrid = min(cminusgrid);
cminusgrid = cminusgrid + abs(mincminusgrid);
cminusgrid = sign(cminus).*cminusgrid;
[cminusgrid,Iunique] = unique(cminusgrid);
cminus = cminus(Iunique);

cpercent = -200:40:200;
cgridpercent = cpercent;

vmrzon1 = log(vmrzon1);
vmrzon2 = log(vmrzon2);
% vmrzondif = log(vmrzondif);
signvmrzondif = sign(vmrzondif);
vmrzondif = log(abs(vmrzondif));
vmrzondif = vmrzondif + abs(mincminusgrid);
vmrzondif = signvmrzondif.*vmrzondif;

% vmrzondif = -5*ones(size(vmrzondif));

test = cminusgrid;
test1 = vmrzon1;
test2 = vmrzondif;

ngas = length(vmrzon1(1,1,:));
fs = 12;
ylim2 = 300;
ylim1 = 10^-1;
axx1 = nan(ngas,1);
ytickspace = [10^-4 10^-2 1 10^2];
% figpos = [-1262 41  353  893];
figpos = [-1156 40  300  893];
% figpos = [97  -214   480   893];



if yplot == 1
    % Do for the zonal vmr
    figi = randi(100);
    %     figure(figi), set(gcf,'Position', [5,12,1096,704])
    %     figure(figi), set(gcf,'Position', [358,61,722,532])
    %     figure(figi), set(gcf,'Position', [97,49,852,630])
    figure(figi), set(gcf,'Position', figpos)
    %     suptitle(sprintf('%s zonal VMR', gasname))
    %     hold on
    fignames = {'ClO\_meas', 'HOCl\_meas', 'HCl\_meas', 'ClONO2\_meas', 'Cly\_meas'};
    for i = 1:ngas
        axx1(i) = subplot(ngas,1,i);
        contourf(axx1(i), latace, pace, vmrzon1(:,:,i), cgrid); caxis([cgrid(1),cgrid(end)])
        title(fignames{i})
        set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','ON','XMinorTick','on', 'XMinorGrid', 'ON', 'YGrid', 'ON', 'FontSize',fs);
        yticks(ytickspace); ylim([ylim1 ylim2]);
%         if i ~= 5
%             ylabel('')
%             set(gca,'yticklabel',{[]})
%         end
        if i ~= 5
            xlabel('')
            set(gca,'xticklabel',{[]})
        end
        if i == 5
            colorbar('YTick', cgrid, 'Yticklabel', c, 'Location','southoutside')
            xlabel('latitude bins [deg N]')
            ylabel('pressure [hPa]')
            linkaxes(axx1,'xy')
        end
    end
    % second data
    figure(figi + 1), set(gcf,'Position', figpos)
    fignames = {'ClO\_cmam', 'HOCl\_cmam', 'HCl\_cmam', 'ClONO2\_cmam', 'Cly\_cmam'};
    for i = 1:ngas
        axx1(i) = subplot(ngas,1,i);
        contourf(axx1(i), latace, pace, vmrzon2(:,:,i), cgrid); caxis([cgrid(1),cgrid(end)])
        title(fignames{i})
        set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','ON','XMinorTick','on', 'XMinorGrid', 'ON', 'YGrid', 'ON', 'FontSize',fs);
        yticks(ytickspace); ylim([ylim1 ylim2]);
        if i ~= 5
            ylabel('')
            set(gca,'yticklabel',{[]})
        end
        if i ~= 5
            xlabel('')
            set(gca,'xticklabel',{[]})
        end
        if i == 5
            colorbar('YTick', cgrid, 'Yticklabel', c, 'Location','southoutside')
            xlabel('latitude bins [deg N]')
            ylabel('pressure [hPa]')
            linkaxes(axx1,'xy')
        end
    end
    % difference in data
    fignames = {'ClO\_dif', 'HOCl\_dif', 'HCl\_dif', 'ClONO2\_dif', 'Cly\_dif'};
    figure(figi + 2), set(gcf,'Position', figpos)
    for i = 1:ngas
        axx1(i) = subplot(ngas,1,i);
        contourf(axx1(i), latace, pace, vmrzondif(:,:,i), cminusgrid); caxis([cminusgrid(1),cminusgrid(end)])
        title(fignames{i})
        set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','ON','XMinorTick','on', 'XMinorGrid', 'ON', 'YGrid', 'ON', 'FontSize',fs);
        yticks(ytickspace); ylim([ylim1 ylim2]);
        colormap(lbmap(100,'BlueGray'))
%         colormap(bone)
        if i ~= 5
            ylabel('')
            set(gca,'yticklabel',{[]})
        end
        if i ~= 5
            xlabel('')
            set(gca,'xticklabel',{[]})
        end
        if i == 5
            colorbar('YTick', cminusgrid, 'Yticklabel', cminus, 'Location','southoutside')
            xlabel('latitude bins [deg N]')
            ylabel('pressure [hPa]')
            linkaxes(axx1,'xy')
        end
    end
    % percent difference in data
    figure(figi + 3), set(gcf,'Position', figpos)
    fignames = {'ClO\_%dif', 'HOCl\_%dif', 'HCl\_%dif', 'ClONO2\_%dif', 'Cly\_%dif'};
    for i = 1:ngas
        axx1(i) = subplot(ngas,1,i);
        contourf(axx1(i), latace, pace, vmrzondifp(:,:,i), cgridpercent); caxis([cgridpercent(1),cgridpercent(end)])
        title(fignames{i})
        set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','ON','XMinorTick','on', 'XMinorGrid', 'ON', 'YGrid', 'ON', 'FontSize',fs);
        yticks(ytickspace); ylim([ylim1 ylim2]);
        colormap(lbmap(100,'BlueGray'))
%         colormap(bone)
        if i ~= 5
            ylabel('')
            set(gca,'yticklabel',{[]})
        end
        if i ~= 5
            xlabel('')
            set(gca,'xticklabel',{[]})
        end
        if i == 5
            colorbar('YTick', cgridpercent, 'Yticklabel', cpercent, 'Location','southoutside')
%             colorbar('Location','southoutside')
            xlabel('latitude bins [deg N]')
            ylabel('pressure [hPa]')
            linkaxes(axx1,'xy')
        end
    end
    
end
%
end


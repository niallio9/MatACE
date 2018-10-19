function [ ] = plot_ace_climatology_comparison( vmrzon_dif, vmrzonvar_dif, obscount_dif, gas_in )
%A funcion to compare the climatology made by Jaho, with the current
%version of the climatology. The assumption is that the two versions are
%made on the same latitude and altitude grid, and ar for the same gas.


% *INPUT*
%           vmrzon_dif: CELL OF ARRAYS - the difference of the zonal vmrs
%           for each month of the year. 
%
%           vmrzonvar_dif: CELL OF ARRAYS - the difference of the standard
%           deviations of the zonal vmrs for each month of the year. 
%
%           obscount_dif: CELL OF ARRAYS - the difference of the
%           observation counts for each month of the year. 
%
% *OUTPUT*
%           makes a series of plots 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define some things
fs = 16;
ytickspace = [10^-4 10^-2 1 10^2];
% yticknames = ['10^-4','10^-2',]
lat = -87.5:5:87.5;
plev = [1000 850 700 500 400 300 250 200 170 150 130 115 100 90 80 70 50 30 20 ...
    15 10 7 5 3 2 1.5 1 0.7 0.5 0.3 0.2 0.15 0.1 0.08 0.05 0.03 0.02 0.01 0.007 0.004 ...
    0.003 0.002 0.001 0.0008 0.0005 0.0003 0.0002 0.0001]';
monthnames = {'January', 'February', 'March', 'April', 'May', 'June', ...
    'July', 'August', 'September', 'October', 'November', 'December'};

%% do the plots for each variable
% zonal vmr
figi = randi(100);
figure(figi), set(gcf,'Position', [5,12,1096,704]),
hold on
if nargin > 3
    gasname = gas_in;
    title_out = sprintf('difference in %s zonal VMR', gasname);
else
    title_out = sprintf('difference in zonal VMR');
end
suptitle(title_out,0.98, 20)
cmax = 1e-6;
for i = 1:12 % make the subplots for each month
    subplot(3,4,i), pcolor(lat, plev, vmrzon_dif{i});
    set(gca, 'Ydir','reverse', 'YScale', 'log')%, 'YMinorTick','on');
    yticks(ytickspace);
    yticklabels()
    % make the subplots for each month
    title(monthnames{i}(1:3))
    cmap = lbmap(80,'BrownBlue');
    colormap(cmap);
    caxis([-cmax,cmax]);
%     h1.LevelStep = 0.2;
%     h1.ShowText = 'on';
    set(gca,'FontSize',fs)
end
%
end


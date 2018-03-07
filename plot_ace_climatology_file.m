function [ ] = plot_ace_climatology_file( climfile_in )
%A function to plot the climatology data for ACE. The function uses
%'pcolor'.
%
% *INPUT*
%           climstruct_in: STRING - the name of the file containing the ACE
%           climatology data.  
%
% *OUTPUT*
%           makes a plot 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define some things and plot
load(climfile_in); % loads a variable called climstruct
clim = climstruct;
plot_ace_climatology(clim);
%
end


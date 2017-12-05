function [medi,meani,maxi,mini,vari] = get_LST_stats(LST)
% this code is slightly modified from an original script by Jaho called
% 'make_LST.m'. -NJR 11/17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% program that creates the relevant information for the LST for a lat/alt vector of data
% Input
% LST   = Local Solar Time: a vector
% Output
% medi = median
% meani = mean
% maxi  = max
% mini  = min
% stdi  = standard deviation related to the mean
% mdti  = median deviation 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LST(LST==0)=nan;

medi = nanmedian(LST);
meani = nanmean(LST);
maxi  = max(LST);
mini  = min(LST);
vari  = nanvar(LST);
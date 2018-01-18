function [tmjd, dates] = cmam2mjd(cmam_time)
%% this function converts the date in CMAM data, which is 'days since 1900-01-01', to MJD
%% NJR - 08/2016

tcmam = cmam_time;

tmjd = date2mjd(1900,1,tcmam);

%% change the dates to acccount for all of the previuos leap years
dates = datevec(mjd2datenum(tmjd));
minyear = min(dates(:,1));
addon = floor((minyear - 1900)/4 + 1);
tmjd = tmjd + addon;


%% remove the leap dates and correct the mjd

dates = datevec(mjd2datenum(tmjd)); 

I1 = find(ismember(dates(:,2:3), [2,29], 'rows'));

for i = 1:4:length(I1) % jumps in 4 because there are 4 measurements each day
    I = find(ismember(dates(:,2:3), [2,29], 'rows'));
    tmjd(I(1):end) = tmjd(I(1):end) + 1;  % from the first index, add a day on to the rest of the data
end

dates = datevec(mjd2datenum(tmjd));

end
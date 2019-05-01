function [lst, EOT]=mjd2lst(mjd, longitude, do_EOT)
% modified from Jaho's utc2lst
% NJR - 11/17 
% NJR - 03/18 - include EOT correction

if isvector(mjd) && ismatrix(longitude) % for the code to run on stupid Deluge
    mjd = repmat(mjd, length(longitude(:,1)), 1);
elseif isvector(longitude) && ismatrix(mjd)
    longitude = repmat(longitude, 1, length(mjd(1,:)));
end

[~,~,~,hour,minute,secs,ticks] = mjd2utc(mjd);
% whos
lst = mod(hour+minute/60+(secs+ticks)/60/60+longitude/15.0, 24.0);

%% this is a new part to include the Equation Of Time: https://www.spaceacademy.net.au/library/notes/eot.htm
if nargin < 3
    do_EOT = 1;
end

if do_EOT == 1
%     disp('applying EOT correction to LST...')
    if ismatrix(mjd)
        d = nan(size(mjd));
        for i = 1:length(mjd(:,1))
            d(i,:) = mjd2doy(mjd(i,:));
        end
    else
        d = mjd2doy(mjd);
    end
    
    E1 = -7.64 * sin(pi*(d - 2)/180);
    E2 = 9.86 * sin(2*pi*(d - 80)/180);
    EOT = (E1 + E2)/60;
    
    lst = lst + EOT;
%     disp('done')
else 
    EOT = 0;
end
 



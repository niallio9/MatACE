function [ apriori_out, tanstruct_out ] = scale_apriori( apriori_in, tanstruct_in, altitude_or_pressure, limits )
%A function to scale an a priori vmr profile to approximate the input data
%profiles, and add the scaled data to areas outside the pressure limits of
%the input data.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. The GLC data must also be added to
%           the tanstruct so that it has the latitude information.
%
%           apriori_in: ROW VECTOR - a profile of the input gas that will
%           act as the a priori information. The vector should correspond
%           to the altitudes of the tanstruct data. Interpolate to the
%           correct grid if not.
%
%           pressure_hPa_limits: VECTOR - a vector with the lower and upper
%           limits of the useable data.Unites are hPa. e.g., [100, 1].
%
% *OUTPUT* 
%           apriori_out: 2-D ARRAY - scaled a priori values corresponding
%           to the occultations of the input data.
%
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with scaled a priori information added
%           and edited quality flags. This is output if the
%           "pressure_hPa_limits" is entered.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 08/18

%% define some things
prior = apriori_in;
gas = tanstruct_in;
norbit = length(gas.occultation);
prior_out = nan(size(gas.vmr));
factor_limit = 10;
index_limit = 0; 
failure_low = 0;
failure_high = 0;
for n = 1 : norbit
%     n
    inonan = find(~isnan(gas.vmr(:,n)) & gas.vmr(:,n) ~=0); %the indices of vmr values that are not nans
    if ~isempty(inonan)
        %         whos
        factor_low = abs(gas.vmr(inonan(1),n) ./ prior(inonan(1))); % we won't use negative vlaues for the scaled a priori
%         if sign(factor_low) == -1
%             factor_low = 1;
%         end
        if abs(factor_low) < factor_limit
            prior_out(1:inonan(1)-1, n) = prior(1:inonan(1)-1) * factor_low; % the scaled apriori values for the lower altitudes
        else
            fprintf('\nlow altitude factor is > %f. norbit = %i\n', factor_limit, n)
            failure_low = failure_low + 1;
%             diff_low = gas.vmr(inonan(end),n) - prior(inonan(end))
%             prior_out(1:inonan(1)-1, n) = prior(1:inonan(1)-1) + diff_low;
        end
        if inonan(end) >= index_limit
            factor_high = abs(gas.vmr(inonan(end),n) ./ prior(inonan(end))); % we won't use negative vlaues for the scaled a priori
            if abs(factor_high) < factor_limit
                prior_out(inonan(end)+1:end, n) = prior(inonan(end)+1:end) * factor_high; % the same for the higher altitudes
            else
                fprintf('\nhigh altitude factor is > %0.2f. norbit = %i\n', factor_limit, n)
                failure_high = failure_high + 1;
                %             diff_high = gas.vmr(inonan(end),n) - prior(inonan(end))
                %             prior_out(inonan(end)+1:end, n) = prior(inonan(end)+1:end) + diff_high;
            end
        else
            fprintf('\nhighest index of available data (%i) is lower than index limit of %i. norbit = %i\n', inonan(end), index_limit, n)
        end
    else
        % the prior stays as nans for this occultation
    end
end

%% if there are pressure limits entered
if nargin > 2
    disp('adding scaled a priori data to tanstruct')
    gasout = gas;
    switch altitude_or_pressure
        case 'pressure'
            disp('using entered pressure limits (hPa)')
            % fill in the values outside the limits with the scaled apriori
            limit_low = limits(1);
            limit_high = limits(2);
            for n = 1 : norbit
                %low altitude values
                ilow = find(gasout.pressure_hPa(:,n) > limit_low);
                gasout.vmr(ilow,n) = prior_out(ilow,n);
                gasout.vmr_error(ilow,n) = -888;
                gasout.quality_flags(ilow,n) = 8;
                %high altitude values
                ihigh = find(gasout.pressure_hPa(:,n) < limit_high);
                gasout.vmr(ihigh,n) = prior_out(ihigh,n);
                gasout.vmr_error(ihigh,n) = -888;
                gasout.quality_flags(ihigh,n) = 8;
            end
        case 'altitude'
            disp('using entered altitude limits (km)')
            % fill in the values outside the limits with the scaled apriori
            limit_low = limits(1);
            limit_high = limits(2);
            for n = 1 : norbit
                %low altitude values
                ilow = find(gasout.altitude_km < limit_low);
                gasout.vmr(ilow,n) = prior_out(ilow,n);
                gasout.vmr_error(ilow,n) = -888;
                gasout.quality_flags(ilow,n) = 8;
                %high altitude values
                ihigh = find(gasout.altitude_km > limit_high);
                gasout.vmr(ihigh,n) = prior_out(ihigh,n);
                gasout.vmr_error(ihigh,n) = -888;
                gasout.quality_flags(ihigh,n) = 8;
            end
    end
    tanstruct_out = gasout;
    
end

apriori_out = prior_out;
disp(failure_low)
disp(failure_high)
%
end


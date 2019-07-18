function [ tanstruct_1_out, tanstruct_2_out ] = combine_ace_profiles( tanstruct_1_in, tanstruct_2_in, altitude_or_pressure, limits )
%A function to merge data by occultation. The data from tanstruct_2 is
%added to tanstruct_1 where there is are NaNs in tanstruct_1, or according
%to the provided limits. It is assumed that the data sets are on the same
%altitude or pressure grid.

% *INPUT*
%           tanstruct_1_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. The GLC data must also be added to
%           the tanstruct so that it has the latitude information.
%
%           tanstruct_2_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. The GLC data must also be added to
%           the tanstruct so that it has the latitude information.
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
%   NJR - 05/19

%% user defined
% when no altitude/pressure limits are input, this defines whether to apply
% the function to the NaNs found at alitudes above or below where there is
% available data. Options are 'high', 'low', or 'both'
high_or_low = 'both';


%% define some things
[gas1, gas2] = match_ace_data(tanstruct_1_in, tanstruct_2_in);

% make sure that the altitudes match 
if ~isequal(gas1.altitude_km, gas2.altitude_km) && ~isequal(gas1.pressure_hPa, gas2.pressure_hPa)
    error('either altitude or pressure fields for the input data sets must match match. Stopping...')
end
norbit = length(gas1.occultation);
gasout = gas1;
if nargin > 2
    switch altitude_or_pressure
        case 'pressure'
            disp('using entered pressure limits (hPa)')
            % fill in the values outside the limits with the scaled apriori
            limit_low = limits(1);
            limit_high = limits(2);
            for n = 1 : norbit
                %low altitude values
                ilow = find(gasout.pressure_hPa(:,n) > limit_low);
                gasout.vmr(ilow,n) = gas2.vmr(ilow,n);
                gasout.vmr_error(ilow,n) = gas2.vmr(ilow,n);
                gasout.quality_flags(ilow,n) = gas2.quality_flags(ilow,n);
                %high altitude values
                ihigh = find(gasout.pressure_hPa(:,n) < limit_high);
                gasout.vmr(ihigh,n) = gas2(ihigh,n);
                gasout.vmr_error(ihigh,n) = gas2.vmr_error(ihigh,n);
                gasout.quality_flags(ihigh,n) = gas2.quality_flags(ihigh,n);
            end
        case 'altitude'
            disp('using entered altitude limits (km)')
            % fill in the values outside the limits with the scaled apriori
            limit_low = limits(1);
            limit_high = limits(2);
            for n = 1 : norbit
                %low altitude values
                ilow = find(gasout.altitude_km < limit_low);
                gasout.vmr(ilow,n) = gas2.vmr(ilow,n);
                gasout.vmr_error(ilow,n) = gas2.vmr_error(ilow,n);
                gasout.quality_flags(ilow,n) = gas2.quality_flags(ilow,n);
                %high altitude values
                ihigh = find(gasout.altitude_km > limit_high);
                gasout.vmr(ihigh,n) = gas2.vmr(ihigh,n);
                gasout.vmr_error(ihigh,n) = gas2.vmr_error(ihigh,n);
                gasout.quality_flags(ihigh,n) = gas2.quality_flags(ihigh,n);
            end
    end
    
else
    switch high_or_low
        case 'both'
            disp('using limits defined by NaN data at all altitudes')
        case 'high'
            disp('using limits defined by NaN data at higher altitudes')
        case 'low'
            disp('using limits defined by NaN data at lower altitudes')
    end
    for n = 1 : norbit
        %     n
        inonan = find( ~isnan(gas1.vmr(:, n)) ); %the indices of vmr values that are not nans
        if ~isempty(inonan)
            switch high_or_low
                case 'low'
                    gasout.vmr(1:inonan(1) - 1, n) = gas2.vmr(1 : inonan(1) - 1, n); % the profile values for the lower altitudes
                    gasout.vmr_error(1:inonan(1) - 1, n) = gas2.vmr_error(1 : inonan(1) - 1, n);
                    gasout.quality_flags(1 : inonan(1) - 1, n) = gas2.quality_flags(1 : inonan(1) - 1, n);
                case 'high'
                    gasout.vmr(inonan(end) + 1 : end, n) = gas2.vmr(inonan(end) + 1 : end, n); % the same for the higher altitudes
                    gasout.vmr_error(inonan(end) + 1 : end, n) = gas2.vmr_error(inonan(end) + 1 : end, n);
                    gasout.quality_flags(inonan(end) + 1 : end, n) = gas2.quality_flags(inonan(end) + 1 : end, n);
                case 'both'
                    gasout.vmr(1:inonan(1) - 1, n) = gas2.vmr(1 : inonan(1) - 1, n); % the profile values for the lower altitudes
                    gasout.vmr_error(1:inonan(1) - 1, n) = gas2.vmr_error(1 : inonan(1) - 1, n);
                    gasout.quality_flags(1 : inonan(1) - 1, n) = gas2.quality_flags(1 : inonan(1) - 1, n);
                    gasout.vmr(inonan(end) + 1 : end, n) = gas2.vmr(inonan(end) + 1 : end, n); % the same for the higher altitudes
                    gasout.vmr_error(inonan(end) + 1 : end, n) = gas2.vmr_error(inonan(end) + 1 : end, n);
                    gasout.quality_flags(inonan(end) + 1 : end, n) = gas2.quality_flags(inonan(end) + 1 : end, n);
            end
        end
    end
end

tanstruct_1_out = gasout;
tanstruct_2_out = gas2;
%
end


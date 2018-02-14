function [ tanstruct_out ] = merge_ace_dmp( tanstruct_in, dmpstruct_in)
%A function to read ACE DMP data, according to the occultations of
%the input, and add the latitude and longitutde information to the input
%structure. The input ACE data should use the original ACE 1km-grid. The
%tangent latitude and longitude are used when there is no information in
%the DMP file.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           glcstruct_in: STructure - contains the location data for the
%           ace measurements.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same fields as the
%           input, but with ACE GLC information added: latitude and
%           longitude.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 02/2018
tic
%% Only do if the lat and lon info isnt already there
if ~isfield(tanstruct_in,'lat')
    
    %% Define some things
    gas = tanstruct_in;
    
    gasout = gas;
    dmp = dmpstruct_in;
    
    ldmp = length(dmp.occultation);
    lalt = length(dmp.altitude_km(:,1));
    lorbit = length(gas.occultation);
    lonout = nan(lalt,lorbit);
    latout = nan(lalt,lorbit);
    eqlout = nan(lalt,lorbit);
    


    %% look for the right files according to the info in tanstruct
    fprintf('\nReading the occultation names...\n')
    %% Get the unique codes that identify the occultations/orbit.
    gasorbit(1,1:lorbit) = gas.occultation;
    gasorbit(2,1:lorbit) = gas.sr1ss0;
    dmporbit(1,1:ldmp) = dmp.occultation;
    dmporbit(2,1:ldmp) = dmp.sr1ss0;
    fprintf('Done\n')
    
    %% Check which ones match
    [~,ygas,ydmp] = intersect(gasorbit',dmporbit','rows'); % the indices of where the orbits/occultations match
    
    fprintf('\nReading and adding the DMP data for the given occultations...')
    %fill in the glc data
    for i = 1:length(ygas)
        if sum(ismember(i,ygas)) > 1
            fprintf('orbit number %s', orbitnames(i))
            error('there is more than one match for this occultation number')
        end
        lonout(:,ygas(i)) = dmp.lon(:,ydmp(i));
        latout(:,ygas(i)) = dmp.lat(:,ydmp(i));
        eqlout(:,ygas(i)) = dmp.eql(:,ydmp(i));
    end
    % fill in the rest with the tangent values
    for i = 1:lorbit
        if sum(ismember(i,ygas)) == 0 % if the index does not appear in the intersection (if that orbit is not found in the glc data)
            lonout(:,i) = repmat(gas.lon_tangent(i),lalt,1);
            latout(:,i) = repmat(gas.lat_tangent(i),lalt,1);
            eqlout(:,i) = nan(lalt,1); % just put in nans for the equivalent latitude here. can't just make something up.
        end
    end
    % fill in the nan columns of lat and lon with the tangent values
    for i = 1:lorbit
        if nansum(latout(:,i)) == 0
            lonout(:,i) = repmat(gas.lon_tangent(i),lalt,1);
            latout(:,i) = repmat(gas.lat_tangent(i),lalt,1);
        end
    end
    
    %% add the lat and lon info to the tanstruct
    gasout.lon = lonout;
    gasout.lat = latout;
    gasout.eql = eqlout;
    gasout.eql(gasout.eql == -999) = nan; % change the -999 values to nans
    %% change the size of the other fields to match the altitudes of the DMPs
    % The DMPs follow the ACE 1km altitude grid, but stop at 75.5km (the 76th index of the altitude grid)
    gasout.altitude_km = gasout.altitude_km(1:lalt,:);
    gasout.vmr = gasout.vmr(1:lalt,:);
    gasout.vmr_error = gasout.vmr_error(1:lalt,:);
    gasout.quality_flags = gasout.quality_flags(1:lalt,:);
    gasout.pressure_hPa = gasout.pressure_hPa(1:lalt,:);
    
    tanstruct_out = gasout;
    fprintf('\nDone\n')
else
    tanstruct_out = tanstruct_in;
    fprintf('\nThe lat and lon fields are already in this structure, so nothing was done here.\n');
end
%
toc
end


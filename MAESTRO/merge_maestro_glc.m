function [ tanstruct_out ] = merge_maestro_glc( tanstruct_in, glcstruct_in)
%A function to read ACE GLC data files, according to the occultations of
%the input, and add the latitude and longitutde information to the input
%structure. The input ACE data should use the original ACE 1km-grid. The
%tangent latitude and longitude are used when there is no information in
%the GLC file.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           glcstruct_in: STRUCTURE - contains the location data for the
%           maestro measurements. This is also refered to a a SunriseTable
%           or SunsetTable or SunriseSunsetTable or sr1ss0_table
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same fields as the
%           input, but with MAESTRO GLC information added: latitude and
%           longitude.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017
% NJR - 01/2018
tic
%% Only do if the lat and lon info isnt already there
if ~isfield(tanstruct_in,'lat')
    
    %% Define some things
    gas = tanstruct_in;

    gasout = gas;
    lalt = length(gas.altitude_km(:,1));
    lorbit = length(gas.occultation);
    lonout = nan(lalt,lorbit);
    latout = nan(lalt,lorbit);
    
    glc = glcstruct_in;
    lglc = length(glc.occultation);

    %% look for the right files according to the info in tanstruct
    fprintf('\nReading the occultation names...\n')
    %% Get the unique codes that identify the occultations/orbit.
    gasorbit(1,1:lorbit) = gas.occultation;
    gasorbit(2,1:lorbit) = gas.sr1ss0;
    glcorbit(1,1:lglc) = glc.occultation;
    glcorbit(2,1:lglc) = glc.sr1ss0;
    fprintf('Done\n')
    
    %% Check which ones match
    [~,ygas,yglc] = intersect(gasorbit',glcorbit','rows'); % the indices of where the orbits/occultations match
    
    %The glc files are of the form: ace.ss1234a.txt or ace.ss12345a.txt
    fprintf('\nReading and adding the GLC data for the given occultations...\n')
    %fill in the glc data
    for i = 1:length(ygas)
        if sum(ismember(i,ygas)) > 1
            fprintf('orbit number %s', orbitnames(i))
            error('there is more than one match fot this occultation number')
        end
%         lonout(:,ygas(i)) = glc.lon(:,yglc(i));
%         latout(:,ygas(i)) = glc.lat(:,yglc(i));
% whos
        lonout(:,ygas(i)) = repmat(glc.lon_tangent(yglc(i)),lalt,1);
        latout(:,ygas(i)) = repmat(glc.lat_tangent(yglc(i)),lalt,1);
    end
%     % fill in the rest with the tangent values
%     for i = 1:lorbit
%         if sum(ismember(i,ygas)) == 0 % if the index does not appear in the intersection (if that orbit is not found in the glc data)
%             lonout(:,i) = repmat(gas.lon_tangent(i),lalt,1);
%             latout(:,i) = repmat(gas.lat_tangent(i),lalt,1);
%         end
%     end
%     % fill in the nan columns with the tangent values
%     for i = 1:lorbit
%         if nansum(latout(:,i)) == 0
%             lonout(:,i) = repmat(gas.lon_tangent(i),lalt,1);
%             latout(:,i) = repmat(gas.lat_tangent(i),lalt,1);
%         end
%     end
    
    %% add the lat and lon info to the tanstruct
    gasout.lon = lonout;
    gasout.lat = latout;
    
    tanstruct_out = gasout;
    fprintf('\nDone\n')
else
    tanstruct_out = tanstruct_in;
    fprintf('\nThe lat and lon fields are already in this structure, so nothing was done here.\n');
end
%
toc
end


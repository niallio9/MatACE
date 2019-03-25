function [ tanstruct_out_filtered ] = apply_ace_filter_monthly2X( tanstruct_in)
%A function to create zonally averaged climatologies of ACE measurements,
%by each unique calendar month. 'make_ace_climatology.m' is called here.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. The GLC data must also be added to
%           the tanstruct so that it has the latitude information.
%
% *OUTPUT*
%           tanstruct_out_filtered: STRUCTURE - output has the same
%           fields as the input, but filtered according to the above
%           definition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 02/19

tic
%% Sort some things
data_out = tanstruct_in;
% zlev = length(datain.altitude_km(:,1));
% zlev = length(data_out.pressure_hPa(:,1));
% numel_data_out = numel(data_out);
size_data_out = size(data_out.vmr);

% get the unique years for the data
vdates = datevec(mjd2datenum(data_out.date_mjd));
years_unique = unique(vdates(:,1));
nfound = 0;
%% For each year loop through the months and create climatologies for each.
for j_year = 1:length(years_unique)
%     j_year
    %subset the data by year
    warning off % suppress warnings about reducing the data to zero here. There is output below if this is the case
    datain_yearj = subset_ace_by_year(data_out, years_unique(j_year));
    warning on
    for i_month = 1:12
%         i_month
        %subset the ace data by month
        warning off % suppress warnings about reducing the data to zero here. There is output below if this is the case
        datain_yearj_monthi = subset_ace_by_month(datain_yearj,i_month);
        data_length = length(datain_yearj_monthi.occultation); % number of occultations in the monthly subset
        numel_datain_yearj_monthi = numel(datain_yearj_monthi.vmr); % number of elements in the data
        warning on
        if data_length ~= 0
            %%%% put a while loop here and get the dates of the values that
            %%%% are twice as large as all others for the month. change
            %%%% them to nans
            do_loop = 1;
            while do_loop == 1
%                 disp('loop')
                [xmax, i_xmax] = max(datain_yearj_monthi.vmr(:)); % the linear index of the maximum value
                is2X = find(xmax > 2*datain_yearj_monthi.vmr | isnan(datain_yearj_monthi.vmr));
%                 length(is2X)
%                 size(datain_yearj_monthi.vmr)
%                 numel_datain_yearj_monthi
                if length(is2X) >= numel_datain_yearj_monthi - 1 % if you found a value more than twice as high as all others (except itself, obvs)
                    nfound = nfound + 1;
                    datain_yearj_monthi = remove_ace_data_by_index(datain_yearj_monthi, i_xmax); % remove the value from the month subset
                    % remove the value from the main dataset
                    [I_xmax, J_xmax] = ind2sub(size(datain_yearj_monthi.vmr), i_xmax); % get the row and column indices for the value in the subset
%                     find(data_out.date_mjd ==
%                     datain_yearj_monthi.date_mjd(J_xmax)); 
                    i_xmax_dataout = sub2ind(size_data_out, I_xmax, find(data_out.date_mjd == datain_yearj_monthi.date_mjd(J_xmax))); % the linear index of the point in the main dataset
                    data_out = remove_ace_data_by_index(data_out, i_xmax_dataout);
%                     return
                else
                    do_loop = 0;
                end
            end
        else
            fprintf('There are no ACE data for %i %d\n', i_month, years_unique(j_year));
        end
    end
end
tanstruct_out_filtered = data_out;
fprintf('%i convicts found', nfound)
fprintf('\nDone :)\n')
toc
%
end
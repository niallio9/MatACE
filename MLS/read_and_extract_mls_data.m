function [ mls_out ] = read_and_extract_mls_data( mls_list, swath, save_name )
%A function to read the MLS data in .he5 format and output a cell array of
%structure with the information. Use the 'mls_list' input to subset the files
%to a specific year or gas, etc.

% *INPUT*
%           mls_list: STRING - the name of (or path to) a text file that
%           contains the names of all of the MLS files that you would like
%           to read in. The MLS data generally comes in files that
%           corespond to the day of the year.
%
%           swath: STRING - the name of the gas to which the data
%           corresponds
%
% *OUTPUT*
%           mls_out: CELL ARRAY - each cell corresponds to one .he5 file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 12/18
%%
if nargin < 3
    savename = [];
else
    save_dir = 'C:\Users\ryann\MLS\matdata\'; % this may need to changed sometimes
    savename = save_name;
    if ~isdir(save_dir)
        fprintf('\nIt doesn''t look like ''%s'' exists...\n',save_dir)
        error('The input save directory couldn''t be found')
    end
end
filelist = mls_list;
filein = fopen(filelist);
i = 0;
years = [];
disp('reading filenames...')
while (feof(filein) ~= 1)
    i = i+1;
    filename{i} = fgetl(filein); %don't know how many files are in the list
    yeari = filename{i}(end-11:end-8);
    dayi = filename{i}(end-6:end-4);
    % if the data is a subset
    %       yeari = filename{i}(end-15:end-12);
    %       dayi = filename{i}(end-10:end-8);
    yearandday(i) = str2double(strcat(yeari,dayi)); % gives a unique number of days. the values are sorted
    years = [years,str2double(yeari)];
    years_unique = unique(years);
end
%check if there are any files that correspond to the same year and day.
%Sometimes MLS data has duplicate days for different versions or something,
%idk.
if ~isequal(sort(yearandday),unique(yearandday))
    warning('there are duplicate days in the data')
    yearandday = sort(yearandday);
    diff = yearandday(2:end) - yearandday(1:end-1);
    %     fprintf('the duplicates are for yyyyddd:\n')
    %     yearandday(diff == 0)
    disp('the dates (yyyyddd) of the duplicates have been returned as output')
    mls_out = yearandday(diff == 0)';
    return
end
clear yearandday years
disp('done')

%% set up the structure and read the data from the files
dummy{1,1} = readL2GP(filename{1}, swath);
evalc('dummy = extract_mls_data(dummy);');
if ~strcmp(swath,'ScaledPV')
    lz = length(dummy.vmr(:,1));
else
    lz = length(dummy.spv(:,1));
end
clear dummy
time_est = 1.2e6 * length(years_unique); % there's about 1.1 million files a year
size_est = [lz,time_est];
mls_out.source_file = [];
mls_out.gas = swath;
if ~strcmp(swath,'ScaledPV')
    mls_out.vmr = nan(size_est);
    mls_out.vmr_error = nan(size_est);
else
    mls_out.spv = nan(size_est);
end
mls_out.pressure_hPa = [];
mls_out.date_mjd = nan(1,time_est);
mls_out.lat = nan(1,time_est);
mls_out.lon = nan(1,time_est);


lfiles = length(filename);
fprintf('reading data from the %i files...\n',lfiles)
kend = 0;
for i = 1:lfiles
%     i
    if mod(i,100) == 0
        fprintf('past file %i of %i\n',i,lfiles)
    end
    mls_data{1,1} = readL2GP(filename{i}, swath); %#ok<NASGU>
    evalc('mls_outi = extract_mls_data(mls_data);');
%     mls_outi = extract_mls_data(mls_data);
%     mls_outi
    kstart = kend+1;
    kend = kend + length(mls_outi.date_mjd);
    if ~strcmp(swath,'ScaledPV')
        mls_out.vmr(:,kstart:kend) = mls_outi.vmr;
        mls_out.vmr_error(:,kstart:kend) = mls_outi.vmr_error;
    else
        mls_out.spv(:,kstart:kend) = mls_outi.spv;
    end
    mls_out.date_mjd(:,kstart:kend) = mls_outi.date_mjd;
    mls_out.lat(kstart:kend) = mls_outi.lat;
    mls_out.lon(kstart:kend) = mls_outi.lon;
end
mls_out.pressure_hPa = mls_outi.pressure_hPa;
mls_out.source_file = mls_outi.source_file;
%% sort the data by time
disp('sorting the data by time...')
[mls_out.date_mjd, I] = sort(mls_out.date_mjd);
if ~strcmp(swath,'ScaledPV')
    mls_out.vmr = mls_out.vmr(:,I);
    mls_out.vmr_error = mls_out.vmr_error(:,I);
else
    mls_out.spv = mls_out.spv(:,I);
end
mls_out.lat = mls_out.lat(:,I);
mls_out.lon = mls_out.lon(:,I);
disp('done')

%% clear extra nan columns
disp('removing columns with no data...')
ibad = find(isnan(mls_out.date_mjd));
if ~strcmp(swath,'ScaledPV')
    mls_out.vmr(:,ibad) = [];
    mls_out.vmr_error(:,ibad) = [];
else
    mls_out.spv(:,ibad) = [];
end
mls_out.lat(ibad) = [];
mls_out.lon(ibad) = [];
mls_out.date_mjd(ibad) = [];
disp('done')
mlsstruct = mls_out; %#ok<NASGU>

%% save the data
if ~isempty(savename)
    savedest = fullfile(save_dir,savename);
    fprintf('saving data to %s\n', savedest);
    save(savedest,'mlsstruct','-v7.3');
    fprintf('done\n')
end
%%
disp('all done')
end


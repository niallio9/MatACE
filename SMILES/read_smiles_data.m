function [ smiles_out ] = read_smiles_data( smiles_list, swath )
%A function to read the SMILES data in .he5 format and output a cell array of
%structure with the information. Use the 'mls_list' input to subset the files
%to a specific year or gas, etc.

% *INPUT*    
%           smiles_list: STRING - the name of (or path to) a text file that
%           contains the names of all of the SMILES files that you would like
%           to read in. The SMILES data generally comes in files that
%           corespond to the day of the year.
%
%           swath: STRING - the name of the gas to which the data
%           corresponds
%
% *OUTPUT*
%           smiles_out: CELL ARRAY - each cell corresponds to one .he5 file. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 06/18

filelist = smiles_list;
filein = fopen(filelist);
i = 0;
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
end
%check if there are any files that corespond to the sme year and day.
%Sometimes MLS data has duplicate days for different versions or something,
%idk.
if ~isequal(sort(yearandday),unique(yearandday))
    warning('there are duplicate days in the data')
    yearandday = sort(yearandday);
    diff = yearandday(2:end) - yearandday(1:end-1);
    fprintf('the duplicates are for yyyyddd:\n')
    yearandday(diff == 0)
    return
end
disp('done')
lfiles = length(filename);
smiles_out = cell(1,lfiles);
doy = nan(1,lfiles);
fprintf('reading data from the %i files...\n',lfiles)
for i = 1:lfiles
    if mod(i,100) == 0
        fprintf('past file %i of %i\n',i,lfiles)
    end
    smiles_out{1,i} = readL2GP(filename{i}, swath);
    doy(i) = str2double(smiles_out{1,i}.date(6:8));
end
disp('done')
%
end


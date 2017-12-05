%% A script to move all of the ACE dmp files into a central folder.
%The more general script, 'move_files_up_2levels.m', can replace this.

%-NJR 10/17

%% Define some things
% this folder is the top directory of the DMP files. It contains folders
% that correspond to years
dir_dmp = '/home/niall/dmp3.5/v2.0';
dir_out = '/home/niall/dmp3.5/v2.0_all';

%% go through all of the directories in dir_dmp and move the files to dir_out

dir_years = dir(dir_dmp);
dir_years = {dir_years.name}; % get the name of the directories
dir_years = dir_years(3:end); % remove the parent directories

for i = 1 : length(dir_years)
    % get the names of the daily folders
    dir_days = dir(strcat(dir_dmp,'/',dir_years{i})); 
    dir_days = {dir_days.name};
    dir_days = dir_days(3:end);
    for j = 1 : length(dir_days)
        % get the names of the .asc files
        dmpfiledest = strcat(dir_dmp,'/',dir_years{i},'/',dir_days{j},'/*.asc');
        try
            movefile(dmpfiledest, dir_out)  
        catch
        end
    end

    
end

function [] = move_files_up_1level(top_directory, destination_folder_name)
%%
% A script to moves files up 2 directory levels into a main folder in the
% top directory. The name of the folder is given as input.


%-NJR 10/17

%% Set up some things
dir_top = top_directory; %% this folder is the top directory of the files.
dir_out = strcat(destination_folder_name); % the destination folder is in the same directory as the top directory
if isdir(dir_out)
    
    %% go through all of the directories in dir_top and move the files to dir_out
    
    dir_level1 = dir(dir_top);
    dir_level1 = {dir_level1.name}; % get the name of the directories
    dir_level1 = dir_level1(3:end); % remove the parent directories
    
    
    for j = 1 : length(dir_level1) % loop through the 1st level directories
        % get the names of all the files in this folder
        filedest = dir(strcat(dir_top,'/',dir_level1{j})); % change this line to specify only certain files
        filedest = {filedest.name};
        filedest = filedest(3:end); % remove the parent directories
        for k = 1: length(filedest) % loop throught the files
            %try
            if mod(k,100) == 0
                fprintf('%d of %d\n',k,length(filedest));
            end
            filei = strcat(dir_top,'/',dir_level1{j},'/',filedest{k});
%             copyfile(filei, dir_out)
            movefile(filei, dir_out)
            %catch
            %end
        end
    end
        
        
else
    fprintf('\nthe output directory ''%s'' does not exist. create it first\n',dir_out);
end

end

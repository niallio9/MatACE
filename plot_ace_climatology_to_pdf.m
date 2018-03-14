function [ ] = plot_ace_climatology_to_pdf( climatology_directory )
%A function to plot the climatology data for ACE and save the figures in a
%pdf file for each gas and type of climatology. The function assumes that
%the 'climatology_directory' contains folders for each/any gas. The monthly
%and seasonal climatologies are in the folders specific to each gas, and
%the serial-monthly climatologies are in a subdirectory called
%'serial_month'.
%
% *INPUT*
%           climatology_directory: STRING - the name of the directory
%           containing the ACE climatology data.
%
% *OUTPUT*
%           makes pdf files
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% get the name of the folders in the climatology directory
climdir = climatology_directory;
if ~isdir(climdir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',climdir)
    error('The input directory couldn''t be found')
end
temp = dir(climdir);
temp = {temp.name};
% remove folders that start with '.'
gasfolder = cell(150,1); % will remove the empty cells later
for i = 1:length(temp)
    if ~strcmp(temp{i}(1), '.')
        gasfolder{i} = temp{i};
    end
end
gasfolder = gasfolder(~cellfun('isempty',gasfolder)); % remove the empty cells

%%
for i = 1:length(gasfolder) % loop through the gases
    gasdir_i = fullfile(climdir,gasfolder{i}); % make the full path to the folder for a gas
    templat = dir(fullfile(gasdir_i,'*lat*.mat')); % get the names of the latitude climatology files
    tempeql = dir(fullfile(gasdir_i,'*eql*.mat')); % get the names of the equivalent latitude climatology files
    templat_serial = dir(fullfile(gasdir_i,'serial_month','*lat*.mat')); % get the names of the serial-month latitude climatology files
    tempeql_serial = dir(fullfile(gasdir_i,'serial_month','*eql*.mat')); % get the names of the serial-month latitude climatology files
    climfile_lat = {templat.name};
    climfile_eql = {tempeql.name};
    climfile_lat_serial = {templat_serial.name};
    climfile_eql_serial = {tempeql_serial.name};
    
    % make the pdf for the lat files
    if ~isempty(climfile_lat)
        pdf_name = sprintf('%s/%s_lat_climatologies.pdf', climdir, gasfolder{i}); % the name of the appended pdf (output)
        for j = 1:length(climfile_lat)
            climfile_i = fullfile(gasdir_i, climfile_lat{j}); % the name of a climatology file
            plot_ace_climatology_file(climfile_i);
            pdf_file_i = strcat(climfile_i(1:end-4),'.pdf'); % make pdf files with the same name as the mat files
            export_fig(pdf_file_i);
            close all
        end
        % append the pdf files into one pdf file
        temppdf = dir(fullfile(gasdir_i,'*lat*.pdf'));
        pdf_file_lat = fullfile(gasdir_i,{temppdf.name});
        if exist(pdf_name,'file') == 2 % check if a previous appended pdf file exists and delete it if so
            delete(pdf_name)
        end
        append_pdfs(pdf_name,pdf_file_lat{:});
        delete(fullfile(gasdir_i,'*lat*.pdf')); % delete the temporary pdf files
    end
    
    % make the pdf for the eql files
    if ~isempty(climfile_eql)
        pdf_name = sprintf('%s/%s_eql_climatologies.pdf', climdir, gasfolder{i}); % the name of the appended pdf (output)
        for j = 1:length(climfile_eql)
            climfile_i = fullfile(gasdir_i, climfile_eql{j}); % the name of a climatology file
            plot_ace_climatology_file(climfile_i);
            pdf_file_i = strcat(climfile_i(1:end-4),'.pdf'); % make pdf files with the same name as the mat files
            export_fig(pdf_file_i);
            close all
        end
        % append the pdf files into one pdf file
        temppdf = dir(fullfile(gasdir_i,'*eql*.pdf'));
        pdf_file_eql = fullfile(gasdir_i,{temppdf.name});
        if exist(pdf_name,'file') == 2 % check if a previous appended pdf file exists and delete it if so
            delete(pdf_name)
        end
        append_pdfs(pdf_name,pdf_file_eql{:});
        delete(fullfile(gasdir_i,'*eql*.pdf')); % delete the temporary pdf files
    end
    
    % make the pdf for the serial-month lat files
    if ~isempty(climfile_lat_serial)
        pdf_name = sprintf('%s/%s_lat_serial_climatologies.pdf', climdir, gasfolder{i}); % the name of the appended pdf (output)
        for j = 1:length(climfile_lat_serial)
            climfile_i = fullfile(gasdir_i,'serial_month', climfile_lat_serial{j}); % the name of a climatology file
            plot_ace_climatology_file(climfile_i);
            pdf_file_i = strcat(climfile_i(1:end-4),'.pdf'); % make pdf files with the same name as the mat files
            export_fig(pdf_file_i);
            close all
        end
        % append the pdf files into one pdf file
        temppdf = dir(fullfile(gasdir_i,'serial_month','*lat*.pdf'));
        pdf_file_lat = fullfile(gasdir_i,'serial_month',{temppdf.name});
        if exist(pdf_name,'file') == 2 % check if a previous appended pdf file exists and delete it if so
            delete(pdf_name)
        end
        append_pdfs(pdf_name,pdf_file_lat{:});
        delete(fullfile(gasdir_i,'serial_month','*lat*.pdf')); % delete the temporary pdf files
    end
    
    % make the pdf for the serial-month eql files
    if ~isempty(climfile_eql_serial)
        pdf_name = sprintf('%s/%s_eql_serial_climatologies.pdf', climdir, gasfolder{i}); % the name of the appended pdf (output)
        for j = 1:length(climfile_eql_serial)
            climfile_i = fullfile(gasdir_i,'serial_month', climfile_eql_serial{j}); % the name of a climatology file
            plot_ace_climatology_file(climfile_i);
            pdf_file_i = strcat(climfile_i(1:end-4),'.pdf'); % make pdf files with the same name as the mat files
            export_fig(pdf_file_i);
            close all
        end
        % append the pdf files into one pdf file
        temppdf = dir(fullfile(gasdir_i,'serial_month','*eql*.pdf'));
        pdf_file_lat = fullfile(gasdir_i,'serial_month',{temppdf.name});
        if exist(pdf_name,'file') == 2 % check if a previous appended pdf file exists and delete it if so
            delete(pdf_name)
        end
        append_pdfs(pdf_name,pdf_file_lat{:});
        delete(fullfile(gasdir_i,'serial_month','*eql*.pdf')); % delete the temporary pdf files
    end
    %
end


%
end


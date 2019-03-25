function [ ] = make_list( file_extension, savename )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(savename,'w');
list = dir(file_extension);
list = {list.name}';
for row = 1:length(list)
    fprintf(fid,'%s\n',list{row,1});
end
fclose(fid);
end


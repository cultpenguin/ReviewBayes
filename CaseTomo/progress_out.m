% caseTomo_mul
function progress_out(txt,new_file,file)
if nargin<3
    file = 'progress.txt';
end
if nargin<2
    new_file = 0;
end
currentDateTime = datestr(now, 'yyyy-mm-dd HH:MM:SS');

if new_file == 1
    fileID = fopen(file, 'w');
else
    fileID = fopen(file, 'a');
end
fprintf(fileID, '%s: %s\n', currentDateTime, txt);
fclose(fileID);
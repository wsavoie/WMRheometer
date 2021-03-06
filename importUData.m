function [out] = importfile(fileToRead1)
% %IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 21-Mar-2014 14:54:51

fid = fopen(fileToRead1,'rt');

% Import the file
newData1 = importdata(fileToRead1);
% aaa= textscan(fid,'%s %f', 'delimiter', '\t','collectoutput',true);
aaa= textscan(fid,'%q', 'delimiter', '\t');
p = aaa{end}{end};
z= textscan(p,'%q','delimiter','\t');
z = z{1};
pars = zeros(length(z),1);

for i = 1:2:length(z)
%     eval(strcat(z{i},'=',z{i+1},';'));
    eval(strcat('out((i+1)/2)','=',z{i+1},';'));
end
% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end
fclose('all');

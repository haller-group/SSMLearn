function convertLivescript2Markdown(livescript_name)
% Function that converts a livescript to a markdown file for direct preview
% in github. 
% To execute the function, install the add-on:
% Live Script to Markdown Converter
% Link: https://github.com/roslovets/Live-Script-to-Markdown-Converter

% Convert livescript to md
livescript2markdown(livescript_name,'README.md')

% Read md file 
txt = fileread('README.md');

% Replace alternative text of figures
old = append(pwd,'/README_images/');
txt = replace(txt,old,'');

% Add preview statement
%NL = regexp(txt, '[\n]');
preview_statement = append('This is a preview of the livescript `', livescript_name,'`.',newline,newline);
txt = append(preview_statement,txt);

writelines(txt,'README.md');

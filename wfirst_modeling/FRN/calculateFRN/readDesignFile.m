function kristT = readDesignFile(file_dir, filename)
% Read Krist table

if nargin<=2
    startRow = 2;
    endRow = inf;
end

dotsLoc = strfind(filename,'.');
filetype = filename(dotsLoc(end):end);

switch filetype
    case '.txt'
        
        formatSpec = '%13f%13f%16f%16f%16f%16f%16f%f%[^\n\r]';
        
        % Open the text file.
        fileID = fopen(fullfile(file_dir, filename),'r');
        
        dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', ...
            'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, ...
            'ReturnOnError', false, 'EndOfLine', '\r\n');
        for block=2:length(startRow)
            frewind(fileID);
            dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter',...
                '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, ...
                'ReturnOnError', false, 'EndOfLine', '\r\n');
            for col=1:length(dataArray)
                dataArray{col} = [dataArray{col};dataArrayBlock{col}];
            end
        end
    case '.csv'
        delimiter = ',';
        
        % Format for each line of text:
        % For more information, see the TEXTSCAN documentation.
        formatSpec = '%f%f%f%f%f%f%f%f%[^\n\r]';
        
        % Open the text file.
        fileID = fopen(fullfile(file_dir, filename),'r');
        
        % Read columns of data according to the format.
        % This call is based on the structure of the file used to generate this code. If an
        % error occurs for a different file, try regenerating the code from the Import Tool.
        dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for block=2:length(startRow)
            frewind(fileID);
            dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
            for col=1:length(dataArray)
                dataArray{col} = [dataArray{col};dataArrayBlock{col}];
            end
        end
end


fclose(fileID);
kristT = table(dataArray{1:end-1}, 'VariableNames', {...
    'rlamD','r_arcsec','Intensity','Contrast',...
    'coreThruput','PSFpeak','area_sqarcsec','occTrans'});

return

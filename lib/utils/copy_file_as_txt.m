function copy_file_as_txt(fullPathAndName, outputDirectory)
    
    [filepath, fn, ext] = fileparts(fullPathAndName);
    terminalCommand = sprintf('cp %s %s', fullPathAndName, [outputDirectory filesep fn '.txt']);
    [status, cmdout] = system(terminalCommand, '-echo');

end

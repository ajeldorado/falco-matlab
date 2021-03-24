clear

%% Add Path to Utils
%
addpath(genpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))));

%--------------------------------------------------------------------------
% Define test or tests in form of cell array, this is an input to the
% TestEngine class. User can specify a directory as a test, a particular
% file that includes the directory or location. There is other ways of
% defining tests but for the most part these are the two we will be using.
%--------------------------------------------------------------------------
tests = {[fileparts(fileparts(mfilename('fullpath'))), filesep, 'testarchive']};

%% Example 1
% 
% Runs tests without specifying the output file. This will create a text
% file labeled "TestingResults.txt' with test output results in table
% format.
te=TestEngine(tests);
results1=te.runTests(te.largeSuite);
table(results1)

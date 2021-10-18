clear

%% Add Path to falco-matlab

pathToAdd = fileparts(fileparts(fileparts(mfilename('fullpath')))); % falco-matlab directory;
addpath(genpath(pathToAdd));
rmpath(genpath([pathToAdd filesep '.git']))

%%
%--------------------------------------------------------------------------
% Define test or tests in form of cell array, this is an input to the
% TestEngine class. User can specify a directory as a test, a particular
% file that includes the directory or location. There is other ways of
% defining tests but for the most part these are the two we will be using.
%--------------------------------------------------------------------------
tests = {[fileparts(fileparts(mfilename('fullpath'))), filesep, 'testarchive']};

%% Run tests
fnOut = 'TestingResultsUnit.txt';
te = TestEngine(tests);
results = te.runTests(te.largeSuite, fnOut);
table(results)

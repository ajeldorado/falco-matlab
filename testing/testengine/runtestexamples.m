clear variables
clc

%cd /Users/lmarchen/Documents/MATLAB/Matlab/work/unittesting/unittesting/Tests

%% Add Path to Utils
%
%
%addpath /Users/lmarchen/Documents/MATLAB/Matlab/work/unittesting/unittesting/Tests/utils

%--------------------------------------------------------------------------
% Define test or tests in form of cell array, this is an input to the
% TestEngine class. User can specify a directory as a test, a particular
% file that includes the directory or location. There is other ways of
% defining tests but for the most part these are the two we will be using.
%--------------------------------------------------------------------------
% tests={'/Users/lmarchen/Documents/MATLAB/work/unittesting/Tests/testarchive/classTests/SolverTest.m'...
%      ; '/Users/lmarchen/Documents/MATLAB/work/unittesting/Tests/testarchive/classTests'};
 
tests={'/Users/lmarchen/Documents/MATLAB/Matlab/work/unittesting/unittesting/Tests/testarchive/TestSinc.m'};

%% Example 1
% 
% Runs tests without specifying the output file. This will create a text
% file labeled "TestingResults.txt' with test output results in table
% format.
 te=TestEngine(tests);
 results1=te.runTests(te.largeSuite);

%% Example 2
%
% Runs tests and outputs results to a Matlab "tap" file.
 te=TestEngine(tests);
 results2=te.runTests(te.largeSuite,'TestingResults.tap');

%% Example 3
%
% Runs tests and outputs results to a pdf file. It is a detailed file which
% includes command window output text, summury of test, etc.
te=TestEngine(tests);
results3=te.runTests(te.largeSuite,'TestingResults.pdf');

%% Example 4
% 
% Runs tests and ouputs resuts to an html file.
 te=TestEngine(tests);
 results4=te.runTests(te.largeSuite,'TestingResults.html');


%% Example 5
%
% Runs tests and ouputs results to a docx file.
 te=TestEngine(tests);
 results5=te.runTests(te.largeSuite,'TestingResults.docx');


%% Example 6
%
% Runs tests and ouputs results to a txt file.
 te=TestEngine(tests);
 results6=te.runTests(te.largeSuite,'TestingResults.txt');


%% Example 7: Code Coverage Report
%
% Runs tests and generates a code coverage report to an html file format.
te=TestEngine(tests);
sourceCodedir ='/Users/lmarchen/Documents/MATLAB/Matlab/work/unittesting/unittesting/Tests/utils';
results7=te.genCoverageReport(te.largeSuite,sourceCodedir)


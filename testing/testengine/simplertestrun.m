clear all
clc

%cd /Users/lmarchen/Documents/MATLAB/Matlab/work/unittesting/unittesting/Tests

%% Add Path to Utils
%
%
% addpath /Users/lmarchen/Documents/MATLAB/Matlab/work/unittesting/unittesting/Tests/utils

%--------------------------------------------------------------------------
% Define test or tests in form of cell array, this is an input to the
% TestEngine class. User can specify a directory as a test, a particular
% file that includes the directory or location. There is other ways of
% defining tests but for the most part these are the two we will be using.
%--------------------------------------------------------------------------
tests={'/Users/lmarchen/Documents/falco_source/falco-matlab/testing/testarchive'};

%% Example 1
% 
% Runs tests without specifying the output file. This will create a text
% file labeled "TestingResults.txt' with test output results in table
% format.
 te=TestEngine(tests);
 results1=te.runTests(te.largeSuite);
 table(results1)
 
%% Example 2
%
% Runs tests and outputs results to a Matlab "tap" file.
%  te=TestEngine(tests);
%  results2=te.runTests(te.largeSuite,'TestingResults.tap');
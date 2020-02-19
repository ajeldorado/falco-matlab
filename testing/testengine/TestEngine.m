%% TestEngine:
%
% The TestEngine Class uses the Matlab Testing Frame Work to run the 
% specified tests. The supplied tests can be in the form of files or 
% directories (including subdirectories). It works by creating a largeSuite 
% and then gives the user the option to run tests and report results in
% five different modes. Modes: 1. Tapp mode runs tests in "largeSuite"
% and saves to a ".tap" file which is essentially a text file, 2. PDF mode
% runs tests in largeSuite and generates a detailed report to a
% ".pdf" file, 3. HTML mode runs tests in largeSuite and generates a
% detailed report to an ".html" file, 4. DOCX mode runs tests in
% "largeSuite" and saves detailed results to a ".docx" file, and 5. TXT
% mode runs tests in largeSuite and saves results to a ".txt" file. In
% addition the user can run in Coverage Mode where user provides tests and
% the location of the source code being tested, and generates a Coverage
% Report for code in specified folder and subfolder.
%
% Call Example:
%
% Example 1 (TAPP mode)
% 
% tests = {'test1.m','test2.m',...,'testn.m'};
% writefile = 'myResults.tap';
% te=TestEngine(tests);
% results=te.runTests(te.largeSuite,writefile); 
%
% Example 2 (PDF mode)
%
% tests = {'test1.m','test2.m',...,'testn.m'};
% writefile = 'myResults.pdf';
% te=TestEngine(tests);
% results=te.runTests(te.largeSuite,writefile); 
%
% Example 3 (HTML mode)
%
% tests = {'test1.m','test2.m',...,'testn.m'};
% writefile = 'myResults.html';
% te=TestEngine(tests);
% results=te.runTests(te.largeSuite,writefile); 
%
% Example 4 (DOCX mode)
%
% tests = {'test1.m','test2.m',...,'testn.m'};
% writefile = 'myResults.docx';
% te=TestEngine(tests);
% results=te.runTests(te.largeSuite,writefile); 
%
% Example 5 (TXT mode)
%
% tests = {'test1.m','test2.m',...,'testn.m'};
% writefile = 'myResults.tx';
% te=TestEngine(tests);
% results=te.runTests(te.largeSuite,writefile); 
%
% Or
%
% tests = {'test1.m','test2.m',...,'testn.m'};
% te=TestEngine(tests);
% results=te.runTests(te.largeSuite); 
% (writes to UnitTestingResultsTable.txt)
%
% Example 6 (Coverage Mode)
%
% tests = {'test1.m','test2.m',...,'testn.m'};
% sourceCodedir = {
% te=TestEngine(tests);
% results=genCoverageReport(te.largeSuite,sourceCodedir)
%
classdef TestEngine
    properties
        tests
    end
    properties(Dependent)
        largeSuite
    end
    methods
        function te = TestEngine(tests)
            if nargin > 0
                te.tests = tests;
            end
        end
        function ls=get.largeSuite(obj)
            import matlab.unittest.TestSuite
            ls=[];
            for iter=1:numel(obj.tests)
                if exist(obj.tests{iter},'file')==2
                    suite=TestSuite.fromFile(obj.tests{iter});
                    ls=[ls, suite];
                elseif exist(obj.tests{iter},'dir')==7
                    suite=TestSuite.fromFolder(obj.tests{iter},'IncludingSubfolders',true);
                    ls=[ls, suite];
                end
            end
        end
    end
    methods(Static)
        function results=runTests(largeSuite,writefile)
            if nargin>1
                if contains(writefile,'.tap')
                    import matlab.unittest.TestRunner
                    import matlab.unittest.plugins.TAPPlugin
                    import matlab.unittest.plugins.ToFile
                    
                    runner = TestRunner.withTextOutput;
                    plugin = TAPPlugin.producingOriginalFormat(ToFile(writefile));
                    runner.addPlugin(plugin)
                    results = runner.run(largeSuite);
                elseif contains(writefile,'.pdf')
                    import matlab.unittest.TestRunner;
                    import matlab.unittest.plugins.TestReportPlugin;
                    runner = TestRunner.withNoPlugins;
                    plugin = TestReportPlugin.producingPDF(writefile,'IncludingPassingDiagnostics',true,'IncludingCommandWindowText',true);
                    
                    runner.addPlugin(plugin);
                    results = runner.run(largeSuite);
                elseif contains(writefile,'.html')
                    import matlab.unittest.TestRunner;
                    import matlab.unittest.plugins.TestReportPlugin;
                    runner = TestRunner.withNoPlugins;
                    htmlFolder = 'testResults';
                    plugin = TestReportPlugin.producingHTML(htmlFolder,'MainFile',writefile);
                    runner.addPlugin(plugin);
                    results = runner.run(largeSuite);
                elseif contains(writefile,'.docx')
                    import matlab.unittest.TestRunner;
                    import matlab.unittest.plugins.TestReportPlugin;
                    runner = TestRunner.withTextOutput;
                    plugin = TestReportPlugin.producingDOCX(writefile,'IncludingPassingDiagnostics',true,'IncludingCommandWindowText',true);
                    runner.addPlugin(plugin);
                    results = runner.run(largeSuite);
                else
                    import matlab.unittest.TestRunner
                    runner = TestRunner.withTextOutput;
                    results = runner.run(largeSuite);
                    results_table=table(results);
                    writetable(results_table,writefile)
                end
            else
                import matlab.unittest.TestRunner
                runner = TestRunner.withTextOutput;
                results = runner.run(largeSuite);
                results_table=table(results);
                writetable(results_table,'TestingResults.txt')
            end
        end
    end
    methods(Static)
        function results=genCoverageReport(largeSuite,sourceCodedir)
            import matlab.unittest.TestSuite
            import matlab.unittest.TestRunner
            import matlab.unittest.plugins.CodeCoveragePlugin
            runner = TestRunner.withTextOutput;
            runner.addPlugin(CodeCoveragePlugin.forFolder(sourceCodedir,'IncludingSubfolders',true));
            %runner.addPlugin(CodeCoveragePlugin.forFolder(sourceCodedir));
            results = runner.run(largeSuite);
        end
    end   
end
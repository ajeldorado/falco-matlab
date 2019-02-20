#!/usr/local/bin/python

import io
import shlex
import socket
import os, sys
import logging
import StringIO
import argparse
import subprocess
import numpy as np
import ConfigParser
from string import Template
from time import sleep, time
import multiprocessing as mp
from threading import Timer
#from license import CheckLicenseAvailbility
from multiprocessing.pool import ThreadPool

VENDOR = 'Mathworks'
SOFTWARE = 'MATLAB'
MATLAB_WRAPPER = 'FalcoWrapper.m'
DEFAULT_CONFIG = {'HLC': 'config/falco_config_HLC_default.cfg'}


def Getlogger(logfile):
    """
    Function to configure logs.

    Parameters
    ----------

    Returns
    -------
    logger : logging object
        Pointer to logging object
    """
    # Logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # Create consolehandler
    consoleHandler = logging.StreamHandler()
    consoleHandler.setLevel(logging.INFO)

    # create the file handler
    filehandler = logging.FileHandler(logfile, mode = 'a')
    filehandler.setLevel(logging.INFO)

    # create the logging format
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    filehandler.setFormatter(formatter)
    consoleHandler.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(filehandler)
    logger.addHandler(consoleHandler)

    return logger


def GenerateFileName(wdir, prefix, suffix = '', ext = 'txt', num = 0):
    """Generate output report name.

    Parameters
    -----------
    wdir : str
	Working directory	

    prefix : str
	File name prefix

    ext : str
        File name extension

    num : int
        Unique file number	

    Returns
    -------
    filename : str
	Unique file name
    """
    logger.debug('Generating file name')
   
    filename = os.path.join(wdir, '%s_%s_%d%s.%s' %(prefix, 'Run', num, suffix, ext))

    if not os.path.isfile(filename):
        return filename
    else:
        return GenerateFileName(wdir, prefix, suffix, ext, num+1)


def Worker(args):
    """Parallel worker to initiate a new multiprocessing process
    and run MATLAB. It will be terminated if not completed within
    timeout period.

    Parameters
    ----------
    args : tuple
        Tuple of input arguments

    Returns
    -------
    None
    """
    # Unpack input arguments
    (job_id, coro, cfgfile, param_name, param_value, nthreads, falco_path,
                                       wdir, timeout, logdir, logmsg) = args

    # Start time
    t0 = time() 

    # timeout flag
    timed_out = False

    # Define kill
    kill = lambda process: process.terminate()	

    # Number of iterable parameters
    nparams = len(param_name)

    # Append the parameter to cfgfle
    outfile = GenerateFileName(wdir, '%s_%d' %(coro.upper(), job_id), ext = 'cfg', num = 0)	
    WriteConfig(cfgfile, param_name, param_value, outfile)

    # Unique prefix name
    prefix = coro.upper()
    for i in range(nparams):
        prefix += '_' + param_name[i] + '_' + str(param_value[i])
    fname = GenerateFileName(os.path.join(wdir, 'brief'), prefix, suffix = '_config', ext = 'mat')
    prefix = os.path.split(fname)[1].replace('_config', '').replace('.mat', '')    
    logger.info('Job #%d output Run Label : %s' %(job_id, prefix)) 		

    # Matlab command
    matlab_exe = 'matlab'
    matlab_options = ' -nodisplay -nosplash -nodesktop'
    matlab_script = Template(""" "addpath(genpath('$FALCO_PATH'), '-end'); FalcoWrapper('$CORO','$DEFAULTCFG','$CFGFILE',$NTHREADS,'$WDIR','$PREFIX'); quit;" """)
    #matlab_script = Template(""" "addpath(genpath('$FALCO_PATH'), '-end'); quit;" """)

    # Substitute string template
    d = {'FALCO_PATH': falco_path, 'CORO': coro, 'DEFAULTCFG': DEFAULT_CONFIG[coro.upper()],
	'CFGFILE': outfile, 'NTHREADS': nthreads, 'WDIR': wdir, 'PREFIX': prefix}
    #d = {'FALCO_PATH': falco_path}
    script = matlab_script.substitute(d)

    # Construct matlab command
    cmd = matlab_exe + matlab_options
    cmd += ' -r ' + script

    if logmsg:
        # Log Matlab stdout and stderr to log file for each process
        outlog = os.path.join(wdir, logdir, '%s.log' %prefix)
        with io.open(outlog, 'wb') as writer, io.open(outlog, 'rb', 1) as reader:
	    # Run MATLAB
            process = subprocess.Popen(shlex.split(cmd), bufsize = -1, 
			stdout = writer, stderr = subprocess.STDOUT)
    else:
        # Re-direct STDOUT to DEVNULL
        DEVNULL = open(os.devnull, 'wb')
 
     	# Run MATLAB
        process = subprocess.Popen(shlex.split(cmd), bufsize = -1, 
			stdout = DEVNULL, stderr = subprocess.STDOUT)

    # Define timer
    m_timer = Timer(timeout, kill, (process,))
	
    # Wait for process to end or timeout
    status = ''	
    (stdout, stderr) = (None, None)
    try:
        logger.info('Starting Process #%d' %process.pid)
        m_timer.start()
	
        while process.poll() is None:
            sleep(0.2)

        if process.returncode != 0:
            status = 'ABORT'
            logger.info('Process #%d aborted with returncode %d. Skipping to next job.' %(process.pid, process.returncode))
        else:
            status = 'COMPLETE'
            logger.info('Process #%d finished successfully!' %process.pid)		
    finally:
        m_timer.cancel()

    # Elapsed time to run the job
    elapsed_time = time() - t0
    logger.info('Process #%d elapsed time : %6.2f hours' %(process.pid, elapsed_time/3600.))	

    # Send the result back if process has completed successfully
    if stderr:
        logger.debug('Job #%d stdout message : %s' %(job_id, stdout))
        logger.debug('Job #%d stderr message : %s' %(job_id, stderr))
    else:
	if elapsed_time > timeout:
	    status = 'TIMED OUT'
	    logger.info('Process #%d timed out. Process is terminated.' %process.pid)
	
    return (job_id, process.pid, outfile, prefix, elapsed_time, status)


def WriteConfig(cfgfile, param_name, param_value, outfile):
    """Re-write user config file to include iterable parameter for each
    worker.

    Parameters
    ----------
    param_name : list
        List of parameter names

    param_value : list
        List of parameter values

    Returns
    -------

    """
    logger.info('Writing worker level user config file %s' %outfile)

    # Add a dummy section to read config file
    config = StringIO.StringIO()
    config.write('[dummysection]\n')
    config.write(open(cfgfile).read())
    config.seek(0, os.SEEK_SET)

    cp = ConfigParser.ConfigParser()
    cp.optionxform = str
    cp.readfp(config)

    # Iterate through all the parameters
    nparams = len(param_name)
    for i in range(nparams):
        cp.set('dummysection', param_name[i], param_value[i])

    # Remove dummysection heading and write to output file
    cfg = StringIO.StringIO()
    cp.write(cfg)

    cfg.seek(0, 1)

    with open(outfile, 'w') as fd:
        fd.write(cfg.buf.strip('[dummysection]\n'))

    return


class FalcoWrapper(object):
    def __init__(self, coro, cfgfile, itercfg, ncpus, nthreads, falco_path,
                   	proper_path, wdir, timeout = 28800, comments = '', 
		        checklic = False, logdir = 'logs', logmsg = False):
        """Run MATLAB FALCO on Unix/Linux Servers in Parallel.

        Parameters
        ----------
        coro : str
            Coronagrpah type. E.g. - HLC

        cfgfile : str
            User configuration file

        itercfg : str
            Iterable parameters config file

        ncpus : int
            Number of parallel Matlab sessions to start

        nthreads : int
            How many parallel thread.s each MATLAB session can use

        falco_path : str
            Matlab FALCO path

        proper_path : str
            Matlab PROPER path

        wdir : str
            Working directory where all the results will be written

        Returns
        -------
        None
        """
        # Input parameters
        self.coro = coro
        self.cfgfile = cfgfile
        self.itercfg = itercfg
        self.ncpus = int(ncpus)
        self.nthreads = int(nthreads)
        self.falco_path= falco_path
        self.proper_path = proper_path
        self.wdir = wdir
        self.timeout = int(timeout)
	self.comments = comments
	self.checklic = checklic
	self.logdir = logdir
	self.logmsg = logmsg
	self.report_prefix = 'falco_report'

        # Logger
        self.logger = logging.getLogger(__name__)

        return

    def ReadConfig(self, section = 'Survey'):
        """Read input iterable parameters config file.

        Parameters
        ----------
        section : str
	    Config section name	

        Returns
        -------
        None
        """	
        self.logger.info('Reading user configuration file')

        config = ConfigParser.ConfigParser()
        config.optionxform = str
        config.read(self.itercfg)

        itercfg = {}
        if config.has_section(section):
            # Read all items in the section
            for (key, value) in config.items(section):
                (start, step, end) = [float(val) for val in value.split(':')]
		n_steps = (end+step - start) / step
                itercfg[key] = np.linspace(start, end, n_steps)
        else:
            self.logger.error('Iterable config file does not have Survey section. Stopping.')
            raise ValueError('Iterable config file does not have Survey section. Stopping.')

        self.logger.debug(itercfg)

        return itercfg

    def AllComb(self, *input):
        """Generate all combinations of arg1 and arg2. It is
        cartesian product of parameter1 and parameter2.

        Parameters
        ----------
        input : tuple
            Parameter 1 and 2

        Returns
        -------
            : numpy ndarray
        """
        self.logger.info('Generate all combination of iterable parameters')

	# Number of arguments
	nargs = len(input)

        return np.array(np.meshgrid(*input)).T.reshape(-1,nargs)

    def WriteReportHeader(self, outfile):
        """Write output report to summarize run.

        Parameters
        ----------
        None

        Returns
        -------

        """
        self.logger.info('Writing FALCO report header')

        with open(outfile, 'w') as fd:
            fd.write('#################################################################################\n')
            fd.write('## Coronagraph Type : %s\n' %self.coro.upper())
            fd.write('## Iterable parameter file : %s\n' %self.itercfg)
            fd.write('## Number of parallel workers : %d\n' %self.ncpus)
	    fd.write('## Number of threads per worker : %d\n' %self.nthreads)
            fd.write('## Working Directory : %s\n' %self.wdir)
            fd.write('## Notes : %s\n' %self.comments) 
            fd.write('#################################################################################\n')
            fd.write('Job ID\tProcess ID\tConfig File\tPrefix\tElapsed Time\tStatus\n')
	    fd.write('      \t          \t           \t      \t  (hours)   \t         \n')	
            fd.write('#################################################################################\n\n')
        return


    def __call__(self):
        """Run FALCO wrapper.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # Read iter cfg file
        itercfg = self.ReadConfig()

        # Use AllComb function
        inparam = {'NAME' : [], 'VALUE' : []}
	for key,value in itercfg.iteritems():
	    inparam['NAME'].append(key)
            inparam['VALUE'].append(value)
	
	inparam['VALUE'] = self.AllComb(*inparam['VALUE'])

        # Create config and ws directory in the working directory
        config_dir = os.path.join(self.wdir, 'brief')
        ws_dir = os.path.join(self.wdir, 'ws')
        if not os.path.isdir(config_dir):
            os.makedirs(config_dir)

        if not os.path.isdir(ws_dir):
            os.makedirs(ws_dir)

        # Number of jobs to work on?
        n_jobs = inparam['VALUE'].shape[0]
        self.logger.info('Number of jobs to run : %d' %n_jobs)

        # Use multiprocessing ThreadPool to spawn jobs
        if n_jobs < self.ncpus:
            self.ncpus = n_jobs
	    n_cpu_cores = mp.cpu_count()/2
	    self.nthreads = int(n_cpu_cores/self.ncpus)	

        data = [(i, self.coro, self.cfgfile, inparam['NAME'], inparam['VALUE'][i].tolist(),
 	          self.nthreads, self.falco_path, self.wdir, self.timeout, 
				self.logdir, self.logmsg) for i in range(n_jobs)]
      
	# Generate output file name
	report = GenerateFileName(self.wdir, self.report_prefix, ext = 'txt', num = 0)
	self.logger.info('Output report file name : %s' %report)

        # Check if Matlab license is available (JPL floating license only)
	if self.checklic:
	    if 'jpl.nasa.gov' in  socket.getfqdn():
	        self.logger.info('Checking if MATLAB floating license is available (JPL network only)')
		if CheckLicenseAvailbility(VENDOR, SOFTWARE):
		    self.logger.info('Matlab floating license is available. Continuing with Falco run.')
            	else:
                    self.logger.error('No Matlab floating license available. Stopping')
                    return  	

	# Write output file header
        self.WriteReportHeader(report)
 
        self.logger.info('Starting thread pool with %d workers' %self.ncpus) 
	tp = ThreadPool(self.ncpus)
	results = tp.imap_unordered(Worker, data)
	for result in results:
	    with open(report, 'a') as fd:
		(job_id, process_id, outfile, prefix, elapsed_time, timed_out) = result
		fd.write('%d\t%d\t%s\t%s\t%6.2f\t\t%s\n' %(job_id, process_id,
			   os.path.split(outfile)[1], prefix, elapsed_time/3600., timed_out))

	# Terminate the threadpool
	tp.terminate()	

        return

def ValidateInput(**kwargs):
    """Validate input files and directory to check if they exist.

    Parameters
    ----------
    files : str
	Comma seperated list of files

    dirs : str
	Comma seperated list of directories 	
   
    Returns
    -------
    None	
    """
    if kwargs.has_key('files'):
        for fname in kwargs['files'].split(','):
	    if not os.path.isfile(fname):
	        msg = 'Input file %s does not exist. Stopping.' %fname
        	logger.error(msg)
	        raise ValueError(msg) 			
    
    if kwargs.has_key('dirs'):
        for dirname in kwargs['dirs'].split(','):
	    if not os.path.isdir(dirname):
	        msg = 'Input directory %s does not exist. Stopping.' %dirname
        	logger.error(msg)
	        raise ValueError(msg) 			

    return		

def main(coro, cfgfile, iterfile, wdir, comments, ncpus, falco_path,
        proper_path, timeout, script_path, checklic, logdir, logmsg):
    """Main entry point for the routine.

    Parameters
    ----------
    coro : str
        Coronagraph type

    cfgfile : str
        User configuration file

    iterfile : str
        Iterable parameter file

    wdir : str
        Working directory

    comments : str
        User comments for the run

    ncpus : int
        Number of parallel workers

     falco_path : str
        MATLAB FALCO path

    proper_path : str
        MATLAB PROPER path

    timeout : int
        Job timeout in seconds

    script_path : str
	Full path to the FalcoWrapper.py script	

    checklic : bool
	Check if MATLAB license is available? Default is False
	
    logdir : str
	Log directory name. Logs are written to this directory under the working directory		

    logmsg : bool
	Log Matlab output messages to files?

    Returns
    -------
    None
    """
    # Start time
    t0 = time()

    # Get full path of the file
    cfgfile = os.path.abspath(os.path.expanduser(cfgfile))		
    iterfile = os.path.abspath(os.path.expanduser(iterfile))
    wdir = os.path.abspath(os.path.expanduser(wdir))
    falco_path = os.path.abspath(os.path.expanduser(falco_path))
    proper_path = os.path.abspath(os.path.expanduser(proper_path))
    matlab_script = os.path.join(script_path, MATLAB_WRAPPER)

    # Matlab can only use maximum number of cores on the machine and not threads
    # Based on how many ncpus user asked for, have nthreads so that total 
    # does not exceed number of total processor cores. Otherwise MATLAB will
    # abort with Java memory errors
    #n_cpu_cores = mp.cpu_count()/2    # valid for Intel processors with hyperthreading (2 threads per core)
    n_cpu_cores = mp.cpu_count()
    nthreads = int(n_cpu_cores/ncpus)
	
    # Display info
    logger.info('\n')
    logger.info('FALCO Wrapper to run parallel jobs')
    logger.info(' Coronagraph Type : %s' %coro.upper())
    logger.info(' User configuration file : %s' %cfgfile)
    logger.info(' Iterable parameter file : %s' %iterfile)
    logger.info(' Number of parallel workers : %d' %ncpus)
    logger.info(' Number of threads per worker : %d' %nthreads) 	
    logger.info(' Working Directory : %s' %wdir)
    logger.info(' User comments : %s\n' %comments)

    # Validate input files and directories
    ValidateInput(files = '%s,%s,%s' %(matlab_script, cfgfile, iterfile), 
		  dirs = '%s,%s' %(falco_path, proper_path))	

    fw = FalcoWrapper(coro, cfgfile, iterfile, ncpus, nthreads, falco_path,
                proper_path, wdir, timeout, comments, checklic, logdir, logmsg)
    fw()

    logger.info('Total Elapsed Time : %8.2f hours' %((time() - t0)/3600.))

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Python FALCO Wrapper')
    parser.add_argument('coro', type = str, action = 'store', choices = ['hlc'],
                                help = 'Coronagraph Type')
    parser.add_argument('cfgfile', type = str, action = 'store',
								help = 'User configuration file')
    parser.add_argument('iterfile', type = str, action = 'store',
                                help = 'Iterable parameters config file')
    parser.add_argument('wdir', type = str, action = 'store',
                                help = 'Working directory')

    parser.add_argument('--comments', type = str, action = 'store', default = '',
								help = 'Comments for the run')
    parser.add_argument('--ncpus', type = int, action = 'store', default = 6,
								help = 'Number of parallel workers to spawn')
    parser.add_argument('--falco-path', type = str, action = 'store',
                                default = '/Users/ajriggs/Repos/falco-matlab',
                                help = 'Path to Matlab FALCO')
    parser.add_argument('--proper-path', type = str, action = 'store',
                                default = '/Users/ajriggs/Documents/MATLAB/PROPER',
                                help = 'Path to Matlab PROPER')
    parser.add_argument('--timeout', type = int, action = 'store', default = 18000,
                                help = 'Job timeout in seconds')
    parser.add_argument('--checklic', action = 'store_true', default = False,
                                help = 'Check MATLAB license?')
    parser.add_argument('--log-level', type = str, action = 'store', default = 'info',
                                choices = ['debug', 'info', 'warn', 'error'],
                                help = 'Job timeout in seconds')
    parser.add_argument('--logfile', type = str, action = 'store', default = 'falco.log',
                                help = 'Falco wrapper log file name')
    parser.add_argument('--logdir', type = str, action = 'store', default = 'logs',
                                help = 'Log directory name')
    parser.add_argument('--logmsg', action = 'store_true', default = False,
                                help = 'Log Matlab output messages to file?')


    # Parse input parameters
    args = parser.parse_args()
 	
    # Get absolute path to the python script	
    script_path = os.path.split(os.path.abspath(os.path.expanduser(sys.argv[0])))[0]

    # Create working directory if it does not exist
    if not os.path.isdir(args.wdir):
        os.makedirs(os.path.join(args.wdir, args.logdir))

    # Setup logger
    logger = Getlogger(os.path.join(args.wdir, args.logdir, args.logfile))

    if args.log_level == 'debug':
        logger.setLevel(logging.DEBUG)
    elif args.log_level == 'info':
        logger.setLevel(logging.INFO)
    elif args.log_level == 'warn':
        logger.setLevel(logging.WARN)
    else:
        logger.setLevel(logging.ERROR)

    main(args.coro, args.cfgfile, args.iterfile, args.wdir, args.comments,
        args.ncpus, args.falco_path, args.proper_path,
        args.timeout, script_path, args.checklic, args.logdir, args.logmsg)

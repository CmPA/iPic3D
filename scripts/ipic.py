#!/usr/bin/env python

import os
import sys
import subprocess
import socket # gethostname()
import re # regular expression
#from optparse import OptionParser
import getopt
# http://docs.python.org/2/library/collections.html#collections.deque
from collections import deque # double-ended queue
import inspect
#
# useful documentation:
#
# http://effbot.org/zone/python-list.htm
# http://pymotw.com/2/subprocess/
# http://stackoverflow.com/questions/3777301/how-to-call-a-shell-script-from-python-code

def getdims(inputfile):
    # extract dimensions from intput file
    dims = [1, 1, 1]
    pattern = re.compile(r'^\s*([\w]+)\s*=\s*([\w]+)')
    f = open(inputfile)
    for line in f:
        # key, value = line.split('=')
        #pattern.findall(line)
        match = re.search(pattern, line)
        if match:
            var = match.group(1)
            val = match.group(2)
            if var == 'XLEN':
                dims[0]=int(val)
            elif var == 'YLEN':
                dims[1] = int(val)
            elif var == 'ZLEN':
                dims[2] = int(val)
    return dims
    f.close()

def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno

def issue_command(command):
  if(show):
    print ' '.join(command)
  else:
    print '+', ' '.join(command)
    subprocess.call(command);

def issue_shell_command(command):
  if(show):
    print command
  else:
    print '+', command
    os.system(command)

def construct_run_command(args):

    # convert from deque to list for getopts
    args = list(args)

    # set default values
    num_max_threads = 1
    output = 'data'
    inputfile = 'src/inputfiles/GEM.inp'
    hostname = ''
    mpirun = 'mpiexec'
    global system
    if system == 'xeon':
        mpirun = 'mpiexec.hydra' # is this line needed?
        num_max_threads = 4
    elif system == 'mic':
        mpirun = 'mpiexec.hydra'
        # this should be user configurable
        num_max_threads = 50
        hostname = socket.gethostname()
        micnum = 0
        hostname = hostname + '-mic' + str(micnum)

    try:
      opts, args = getopt.getopt(args, 'i:o:s:t:h:', \
        ['input=', 'output=', 'system=', 'threads=', 'host='])
    except getopt.GetoptError, e:
      if e.opt == 'h' and 'requires argument' in e.msg:
        print 'ERROR: -h requires input filename'
      elif e.opt == 'i' and 'requires argument' in e.msg:
        print 'ERROR: -i requires input filename'
      elif e.opt == 'o' and 'requires argument' in e.msg:
        print 'ERROR: -o requires directory name'
      elif e.opt == 't' and 'requires argument' in e.msg:
        print 'ERROR: -t requires max number of threads'
      elif e.opt == 's' and 'requires argument' in e.msg:
        print 'ERROR: -s requires system name (e.g. "mic" or "xeon")'
      else:
        usage()
        sys.exit(-1)

    for o, a in opts:
        if o in ("-h", "--host"):
          hostname = a
        elif o in ("-i", "--input"):
          inputfile = a
        elif o in ("-o", "--output"):
          output = a
          print 'ERROR: -o is not yet supported'
          sys.exit(1)
        elif o in ("-t", "--threads"):
          num_max_threads = int(a)
        elif o in ("-s", "--system"):
          system = a
        #else:
        #  assert False, "unhandled option"

    if len(args)!=0:
      usage();

    # determine num_procs
    dims = getdims(inputfile)
    XLEN = dims[0]
    YLEN = dims[1]
    ZLEN = dims[2]
    num_procs = XLEN*YLEN*ZLEN
    # num_procs = 4

    arguments = ['exec/iPic3D', inputfile];
    options = ['-n', str(num_procs)]
    if hostname!="":
        options.extend(['-host', hostname])

    if num_max_threads > 1:
        omp_string = 'OMP_NUM_THREADS=' + str(num_max_threads)
        omp = ['-env', omp_string]
        options.extend(omp)

    command = [mpirun]
    command.extend(options)
    command.extend(arguments)
    return command

def ipic_run(args):
    command = construct_run_command(args);
    issue_command(command)

def ipic_show_run(args):
    command = construct_run_command(args);
    print ' '.join(command);

def ipic_make_data():
    # create data subdirectory
    create_data_command = '''mkdir -p data''';
    issue_shell_command(create_data_command)

def ipic_cmake(args):

    # make src a link to the code
    numargs = len(args)
    if numargs==0:
      sourcedir = '..'
    elif numargs==1:
      sourcedir = deque.popleft(args)
    else:
      usage()
      sys.exit()

    if sourcedir!='src':
      rm_command = ['rm -f', 'src'];
      issue_command(rm_command);
      ln_command = ['ln', '-s', str(sourcedir), 'src'];
      issue_command(ln_command)

    ipic_make_data();
    # invoke cmake 
    cmake_command = ['cmake'];
    if system == 'general':
      0
    elif system == 'mic':
      cmake_command.extend(['-DCMAKE_TOOLCHAIN_FILE=src/cmake/cmake_template.cmake.XeonPhi'])
    else:
        print "--system", system, "is not supported"
        sys.exit(-1)
    # issue the command
    cmake_command.extend(['src'])
    issue_command(cmake_command)

def ipic_ctags(args):
    # create tags file using ctags
    create_tags_command = \
        '''find . -name '*.cpp' -or -name '*.h' | grep -v unused | xargs ctags --extra=+qf'''
    issue_shell_command(create_tags_command)
    # sort tags file
    sort_tags_command = '''LC_ALL=C sort -u tags -o tags'''
    issue_shell_command(sort_tags_command)

def ipic_show(args):
    if len(args) == 0:
      ipic_help_show(args)
      sys.exit()
    
    command = deque.popleft(args)
    if command == "run":
      ipic_show_run(args)
    #elif command == "cmake":
    #  ipic_show_cmake(args)
    #elif command == "ctags":
    #  ipic_show_ctags(args)
    else:
        print "ipic show", command, "is not supported"
        sys.exit(-1)

def ipic_basic_help():
    print '''
  To build, you can use:
  
    mkdir build
    cd build
   ''', progname, '''cmake /path/to/ipic3d
    make # or "make -j" to compile in parallel
  
  Then to run the code, use:
  
    ipic run

  If you prefer, use e.g. "ipic show run" to see the shell commands
  that will be executed and then execute them directly yourself.
  
  Available subcommands:

    ''', progname, '''help show    # show what a command would do
    ''', progname, '''help run     # execute iPic3D
    ''', progname, '''help cmake   # execute cmake and create subdirectories
    ''', progname, '''help ctags   # create ctags file to navigate code
    ''', progname, '''help mic     # help for running on mic
    ''', progname, '''help deep    # help for running on deep
  '''

def ipic_help_show(args):
    print '''
  ''', progname, '''show [command]

    show the shell command that would be executed by
      ipic [command]
    '''

def ipic_help_run(args):
    print '''
 ''', progname, '''[-s <mic|xeon>] run [options]

    run iPic3D with appropriate arguments.

    options:
    -t <num_max_threads>: set maximum number of threads
       (default is 1 unless -s <mic|xeon> is set)
    -i <inputfile>: set input file (default is "src/inputfiles/GEM.inp")
    -o <outputdir>: set output directory (default is "data")
    -h <host>: spawn processes on specified host
    '''

def ipic_help_mic(args):
    print '''
  See "ipic help".  Modifications are as follows.

  On the Xeon host processor, use:
  
    ipic -s xeon [command]
  
  On the MIC, use

    ipic -s mic [command]

  To show what a command will do, use e.g.:

    ipic show -s mic [command]
  
  See also:
    ''', progname, '''help deep
    '''

def ipic_help_deep(args):
    print '''
  DEEP needs the following modules:

    module load hdf5/1.8.10-patch1
    module load knc/intel_mpi/4.1.0.030
    module load knc/mic

  For instructions on how to build and run, see
    ''', progname, '''help mic
    '''

def ipic_help_cmake(args):
    print '''
  ''', progname, '''[-s mic] cmake [sourcedir]

     [sourcedir]: the source code directory; by default ".."
     [-s mic]: cross-compile for the mic system
  '''

def ipic_help_ctags(args):
    print '''
  Make sure that you are in the source code directory
  and then run

    ''', progname, '''ctags
    '''

def ipic_help_git(args):
    print '''
    ### This stub gives examples of git commands ###

    # show branch information
    git branch -avv
    # examining the .git directory reveals a wealth of information, e.g.:
    cat .git/config
    # with --stat all files checked in are displayed.
    git log --stat
    # for the following I just do "git tree" (see .gitconfig below):
    git log --oneline --decorate --graph --branches --source
    git status # shows file statuses
    git remote -v # show remote repositories
    # show commits in chronological order.
    git reflog
    # git reflog is useful to get the sha-1 hash of a commit
    # that you recently made and whose branch you accidentally
    # deleted, making it no longer reachable.  Note that
    # each snapshot that you commit should stay in its local
    # repository for 90 days before being garbage collected
    # unless you do something like "git gc".  See also
    # http://gitready.com/advanced/2009/01/17/restoring-lost-commits.html
    #
    # show file
    git show mybranch:myfile
    eg cat myfile # slightly nicer than git show
    # show who checked in what line when under what commit.
    git blame myfile # on current branch
    git blame amaya-library iPic3D.cpp

  for modification:

    # initialize a repository
    mkdir localrepository; cd localrepository; git init
    # creating/removing remote:
    git remote add myremote  https://github.com/alecjohnson/iPic3D.git
    git remote rm myremote  
    # get all branches and their filesystem snapshots
    # from myremote that are not already in localrepository
    git fetch myremote 
    # check in mods
    git stage myfile
    git rm oldfile
    git commit
    # modify a commit message
    git commit --amend
    # create a branch and check it out
    git checkout -b newbranch
    # push branch to server
    eg push --branch newbranch myremote
    # pull changes from server into current branch
    git pull
    # delete a branch on server (!):
    git push myremote --delete mybranch

  # example of global configuration file:

    $ cat ~/.gitconfig
    [user]
    name = eajohnson
    email = e.alec.johnson@gmail.com
    [alias]
    tree = log --oneline --decorate --graph --branches --source
    undo-commit = reset --soft HEAD~1
    '''

def ipic_help(args):
    if len(args) == 0:
      ipic_basic_help()
      sys.exit()
    
    command = deque.popleft(args)
    if command == "show":
      ipic_help_show(args)
    elif command == "run":
      ipic_help_run(args)
    elif command == "mic":
      ipic_help_mic(args)
    elif command == "deep":
      ipic_help_deep(args)
    elif command == "cmake":
      ipic_help_cmake(args)
    elif command == "ctags":
      ipic_help_ctags(args)
    elif command == "git":
      ipic_help_git(args)
    else:
        print "ipic help", command, "is not supported"
        sys.exit(-1)

def usage():
    theline = inspect.currentframe().f_back.f_lineno
    print '  usage() called from ipic.py line ', str(theline)

    print '''
  usage: ''', progname, ''' [show] <command>

  Available commands:
    ''', progname, '''help
    ''', progname, '''show
    ''', progname, '''cmake
    ''', progname, '''ctags
      '''

def ipic_command(argv1):

    global system
    system = 'general'

    # it might be better to use the argparse module rather than getopt,
    # but unfortunately argparse is only available beginning with python 2.7
    # and most HPC platforms seem to have python 2.6 installed.
    # optparse has been deprecated and does not seem to be in python 3;
    # note, however, that argparse was initially an extension of optparse
    # before giving up on backward compatibility.
    #
    try:
      opts, args = getopt.getopt(argv1, 'hs:', ['help', 'system='])
    except getopt.GetoptError, e:
      if e.opt == 's' and 'requires argument' in e.msg:
        print 'ERROR: -s requires system name (e.g. "mic" or "xeon")'
      else:
        usage()
        sys.exit(-1)

    for o, a in opts:
        if o in ("-h", "--help"):
          usage()
          sys.exit()
        elif o in ("-s", "--system"):
          system = a
        #else:
        #  assert False, "unhandled option"

    numargs = len(args)
    if numargs==0:
      usage()
      sys.exit()

    #print args
    args = deque(args)
    command = deque.popleft(args)
    #print list(args)

    if command == "help":
        ipic_help(args)
    # elif command == "show":
    #     ipic_show(args)
    elif command == "ctags":
        ipic_ctags(args)
    elif command == "cmake":
        ipic_cmake(args)
    elif command == "run":
        ipic_run(args)
    else:
        print progname, command, "is not supported"
        sys.exit(-1)

    #print os.path.basename(__file__)
    #print os.path.dirname(__file__)

def main():

    global progname
    progname = os.path.basename(sys.argv[0])
    global dirname
    dirname = os.path.dirname(sys.argv[0])
    global show
    show=0

    argv1 = sys.argv[1:]
    if len(argv1)==0:
      usage()

    if argv1[0]=='show':
      show=1
      argv1=argv1[1:]

    ipic_command(argv1)

if __name__ == '__main__':
    main()


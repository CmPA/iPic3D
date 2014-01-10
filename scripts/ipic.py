#!/usr/bin/env python

import sys
import getopt
# http://docs.python.org/2/library/collections.html#collections.deque
from collections import deque # double-ended queue
import os
#from optparse import OptionParser

# useful documentation:
#
# http://effbot.org/zone/python-list.htm
# http://pymotw.com/2/subprocess/
# http://stackoverflow.com/questions/3777301/how-to-call-a-shell-script-from-python-code

def ipic_ctags(args):
    # create tags file using ctags
    create_tags_command = \
        '''find . -name '*.cpp' -or -name '*.h' | grep -v unused | xargs ctags --extra=+qf'''
    print create_tags_command
    os.system(create_tags_command)
    # sort tags file
    sort_tags_command = '''LC_ALL=C sort -u tags -o tags'''
    print sort_tags_command
    os.system(sort_tags_command)

def ipic_help():
    print '''
  To build, in the iPic3D directory you can use:
  
    rm -rf build # if necessary
    mkdir build
    cd build
    cmake ..
    make # or "make -j" to compile in parallel
  
  To run the code you can use
  
    mkdir data
    mpiexec -n 4 exec/iPic3D ../inputfiles/GEM.inp
  
  where 4 = XLEN times YLEN times ZLEN (defined in GEM.inp).
  
  Available subcommands:

    ''', progname, '''help ctags
    ''', progname, '''help mic
    ''', progname, '''help deep
  '''

def ipic_help_mic(args):
    print '''
  See "ipic help".  Modifications are as follows.

  To run on the Xeon host processor, use something like:
  
    mpiexec.hydra -n 8 -env OMP_NUM_THREADS=4 exec/iPic3D ../inputfiles/GEM.inp
  
  where 8 = XLEN times YLEN times ZLEN.
  
  If you want to cross-compile for the MIC, then the instructions are
  different:
  
      mkdir build.phi
      cd build.phi
      cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/cmake_template.cmake.XeonPhi
      make -j
  
  And to run you use, e.g.:
  
    mkdir data
    mpiexec.hydra -host knc2-mic0 -n 50 -env OMP_NUM_THREADS=4 exec/iPic3D ../inputfiles/GEM.inp
  
  where 50 = XLEN times YLEN times ZLEN.

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

def help(args):
    if len(args) == 0:
      ipic_help()
      sys.exit()
    
    command = deque.popleft(args)
    if command == "mic":
      ipic_help_mic(args)
    elif command == "deep":
      ipic_help_deep(args)
    elif command == "ctags":
      ipic_help_ctags(args)
    elif command == "git":
      ipic_help_git(args)
    else:
        print "ipic help", command, "is not supported"
        sys.exit(-1)

def usage():
    print '''
  usage: ''', progname, ''' [options] <command>

  Available commands:
    ''', progname, '''ctags
    ''', progname, '''help
      '''

def main():

    global progname
    progname = os.path.basename(sys.argv[0])
    global dirname
    dirname = os.path.dirname(sys.argv[0])

    # it might be better to use the argparse module rather than getopt,
    # but unfortunately argparse is only available beginning with python 2.7
    # and most HPC platforms seem to have python 2.6 installed.
    # optparse has been deprecated and does not seem to be in python 3;
    # note, however, that argparse was initially an extension of optparse
    # before giving up on backward compatibility.
    #
    try:
      opts, args = getopt.getopt(sys.argv[1:], 'ho:', ['help', 'output='])
    except getopt.GetoptError, e:
      if e.opt == 'o' and 'requires argument' in e.msg:
        print 'ERROR: -o requires filename'
      else:
        usage()
        sys.exit(-1)

    for o, a in opts:
        if o in ("-h", "--help"):
          usage()
          sys.exit()
        elif o in ("-o", "--output"):
          output = a
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
        help(args)
    elif command == "ctags":
        ipic_ctags(args)
        #print "ctags not yet implemented"
    else:
        print progname, command, "not supported"
        sys.exit(-1)

    #print os.path.basename(__file__)
    #print os.path.dirname(__file__)

if __name__ == '__main__':
    main()


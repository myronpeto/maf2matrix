#!/usr/bin/env python

from __future__ import print_function

import argparse
import os
import subprocess
import sys
from subprocess import *


def execute(cmd):
    print("RUNNING...\n", cmd, "\n")
    process = subprocess.Popen(cmd,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()

    stderr = process.communicate()[1]

    if process.returncode != 0:
        print(
            "[ERROR] command:", cmd, "exited with code:", process.returncode,
            file=sys.stderr
        )
        print(stderr, file=sys.stderr)
        raise RuntimeError
    else:
        return process.returncode


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-m", "--maf-file",
                        type=str,
                        required=True,
                        help="Name of the maf file. Mandatory.")
    parser.add_argument("-d", "--directory",
                        type=str,
                        default="./",
			required=True,
                        help="Working directory. Mandatory.")

    args = parser.parse_args()
    print("maf file: ", args.maf_file, "\ndirectory: ", args.directory)
    
    # JAVA=os.environ["JAVA_HOME"]+"/bin/java"
    
    def jarWrapper(*command):
    	# process = Popen([JAVA, '-jar']+list(command), stdout=PIPE, stderr=PIPE)
    	process = Popen(['java', '-jar']+list(command), stdout=PIPE, stderr=PIPE)
	print('java', '-jar', command[0], command[1], command[2])
	ret = []
    	while process.poll() is None:
        	line = process.stdout.readline()
        	if line != '' and line.endswith('\n'):
            		ret.append(line[:-1])
    	stdout, stderr = process.communicate()
	print(stdout)
    	ret += stdout.split('\n')
    	if stderr != '':
        	ret += stderr.split('\n')
    	ret.remove('')
    	return ret

    command = ['/home/maf2matrix.jar', args.maf_file, args.directory] # Any number of args to be passed to the jar file

    result = jarWrapper(*command)

    print (result)

#!/usr/bin/env python3
"""Compile the primary geometric producer and separate consumer for every profile.

GNU C++17/OpenMP and a little-endian GNU/Linux host are required. Each component
has a 3600 second guard; failed work and logs are retained. Components run
sequentially, with at most four OpenMP threads. Raw witness files are temporary.
"""
from pathlib import Path
import argparse,shutil,sys
from verify_closures import verify,need
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);parser.add_argument('--threads',type=int,default=4,choices=(1,2,3,4));args=parser.parse_args()
 need(shutil.which('g++') is not None and sys.byteorder=='little','GNU C++17 and little-endian host')
 verify(args.output,full=True,threads=args.threads)
if __name__=='__main__':main()

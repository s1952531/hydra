#!/usr/bin/env python3

#=====perform various generic imports=====
import sys,os,shutil,warnings
import numpy as np
#=========================================

#---------------------------------------------------------------------------
# Get list of files for imaging:

filelist=[]
for file in os.listdir("3d"):
   if file.endswith("contours.asc"):
      filelist.append(file)

print ' Choose one of the following files to image:'
print
for i,file in enumerate(filelist):
   print ' ('+str(i+1)+')',file

print
i=int(raw_input(' Choice (default 1): ') or 1)
file=filelist[i-1]
print
print ' Imaging data in',file,'...'
print


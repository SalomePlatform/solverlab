#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2023  CEA
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# 
# See http://www.salome-platform.org or email : webmaster.salome@opencascade.com
# %% LICENSE_END

"""
all file patterns of supposedly (sometimes) created by solverlabGui
replaces '@xxx@' as file.in autotools
"""

"""
https://root.cern.ch/root/html534/src/TRint.cxx.html#OreI4
TRint::ExecLogon()
{
   // Execute logon macro's. There are three levels of logon macros that
   // will be executed: the system logon etc/system.rootlogon.C, the global
   // user logon ~/.rootlogon.C and the local ./.rootlogon.C. For backward
   // compatibility also the logon macro as specified by the Rint.Logon
   // environment setting, by default ./rootlogon.C, will be executed.
   // No logon macros will be executed when the system is started with
   // the -n option
"""

import os
import fnmatch
import getpass
from datetime import datetime
from collections import OrderedDict
import pprint as PP

import xyzpy.loggingXyz as LOG
logger = LOG.getLogger()


#as a dict
_filePatterns = {


###################
"launchSolverlab.sh" : r"""#!/bin/bash

# example of use:
# cd $HOME/SOLVERLABGUI_WORKDIR/solverlabSim_test
# ./launchSolverlab.sh

# this script is in launch directory of solverlab
aDir=`dirname $0`

cd $aDir
echo "launch solverlab.exe code in directory: "`pwd`
which solverlab.exe    #PATH supposed set

if [ -f lock ]; then
    echo "lock file exist and solverlab supposed working in directory."
    echo "delete file lock if you are sure."
    cat lock
    echo "ERROR on locked solverlab directory"
    echo "########## END"
    exit 2
fi

echo "hostname="`hostname` > lock
echo "launchTime="`date` >> lock
echo "launchDirectory="`pwd` >> lock
rm -f lock_old

#launch on file
solverlab.exe etc.

echo "pidSolverlab="$! >> lock
cat lock
#echo "wait"
wait

if [ -s "errors" ]
then
    echo "endTime="`date` >> lock
    mv -f lock lock_old
    echo "ERROR on solverlab"
    cat errors
    echo "solverlab.exe finished"
    echo "########## END"
    exit 2
fi
i

echo "endTime="`date` >> lock
mv -f lock lock_old
echo "solverlab.Exe finished"
echo "########## END"    
""",

"aUserData.dat": r"""#NAME: ExampleName
#TITLE: a data example
#DATE: @DATE@ @TIME@
#COLUMN_NAMES: name1 | name2 | name3 | n__iter__one | n__iter__two
#COLUMN_TITLES: complete name_1 | complete name_2 | complete name_3 | yeah | yooo
#COLUMN_UNITS: cm | cm | cm
#COLUMN_TYPES: D | D | D

1.   2.   3.e+00  0. 4.
1.1  2.1  3.1e+01 1. 3.
1.2  2.2  3.2e+02 2. 2.
1.3  2.3  3.3e+03 3. 1.
1.4  2.4  3.4e+04 4. 0. 

""",


###################
"README.txt" : r'''
###################################
# README.txt @DATE@
# waiting for user modifications...
###################################

...

''',


###################
"aUserFile.sh" : r"""#!/bin/bash

# file @FILE@ as aUserFile.sh
# date @DATE@
# author solverlabGUI @WHOAMI@
# Generated by solverlabGUI, for user use

echo "Begin @FILE@" 

""",


}


###############################################################
def removeDuplicates(*args):
  res = []
  #print "removeDuplicates args",args
  for t in args: res.extend(t)  
  res = list(OrderedDict.fromkeys(res))
  return res
   
###############################################################
def getDefaultDivide(nb):
  """Canvas.Divide by default"""
  if nb <= 1: return (1,1)
  if nb == 2: return (1,2)
  if nb <= 4: return (2,2)
  if nb <= 6: return (2,3)
  if nb <= 9: return (3,3)
  if nb <= 12: return (3,4)
  if nb <= 16: return (4,4)
  nbrow = int(nb/4)
  if nbrow*4 < nb: nbrow += 1
  return (nbrow , 4)

###############################################################
def getTypesOfVisualize():
  check = "visualize_"
  res = [k for k in _filePatterns.keys() if (check in k) and (".C" in k)]
  res = [k.replace(check, "").replace(".C", "") for k in res]
  return sorted(res)

###############################################################
def getVisualizeMethod(aName):
  visualizeTypes = {
    "TDS_graphic": visualize_TDS_graphic,
    "TDS_univariate": visualize_TDS_univariate,
  }
  try:
    return visualizeTypes[aName]
  except:
    logger.error("unknown name: '%s' not in '%s'" % (aName, sorted(visualizeTypes.keys())))
    return None

_tdsNoGarbage = []

###############################################################
def visualize_TDS_graphic(FILEDAT, ATTS, EXPR=[], DIVIDE=None, FILS=[], OPTS=[], Verbose=True):
  """ATTS = ( ("x1:x2"), ("x3:x5") ) for example 2 Draws"""
  import URANIE
  verbose = False
  # Create a TDataServer
  tds = URANIE.DataServer.TDataServer()
  # Read data file
  tds.fileDataRead(FILEDAT)
  if verbose: print("tds %s" % tds)

  #add expressions
  if Verbose: print("****tds add expressions '%s'" % PP.pformat(EXPR))
  for e in EXPR:
    name, formula = e.split("=")
    print("tds.addAttribute '%s'" % e)
    tds.addAttribute(name, formula)

  # Find TCanvas UranieGui
  Canvas = URANIE.ROOT.gROOT.FindObject("UranieGuiCanvas")
  if verbose: print("Canvas %s" % Canvas)
  Canvas.Clear()

  if DIVIDE == None: 
    divi = getDefaultDivide(len(ATTS))
  else:
    divi = DIVIDE

  Canvas.Divide(divi[0], divi[1])
  i = 0
  for att1 in ATTS:
    # Draw filters
    try:
      fil1 = FILS[i]
    except:
      fil1 = ""

    # Draw options
    # see https://root.cern.ch/root/html534/THistPainter.html
    try:
      opt1 = OPTS[i]
    except:
      opt1 = ""

    Canvas.cd(i)
    if Verbose: print("""\
TDS_graphic divide=%s
  attributes='%s'
  filters='%s'
  options='%s'""" % (divi, att1, fil1, opt1))
    tds.Draw(att1, fil1, opt1)
    i += 1
  _tdsNoGarbage.append(tds)
  return tds

###############################################################
def visualize_TDS_univariate(FILEDAT, ATTS, EXPR=[], DIVIDE=None, FILS=[], OPTS=[], Verbose=True):
  import URANIE
  verbose = True
  # Create a TDataServer
  tds = URANIE.DataServer.TDataServer()
  # Read data file
  tds.fileDataRead(FILEDAT)
  if verbose: print("tds %s" % tds)
  
  #add expressions
  print("****tds add expressions '%s'" % PP.pformat(EXPR))
  for e in EXPR:
    name, formula = e.split("=")
    print("tds.addAttribute '%s'" % e)
    tds.addAttribute(name, formula)

  # Find TCanvas UranieGui
  Canvas = URANIE.ROOT.gROOT.FindObject("UranieGuiCanvas")
  if verbose: print("Canvas %s" % Canvas)
  Canvas.Clear()

  if DIVIDE == None: 
    divi = getDefaultDivide(4)
  else:
    divi = DIVIDE
  
  Canvas.Divide(divi[0], divi[1]) 
  # Draw attributes
  att1 = ATTS[0]

  # Draw filters
  fil1 = FILS[0]

  # Draw options
  # see https://root.cern.ch/root/html534/THistPainter.html
  opt1 = OPTS[0]

  if Verbose: print("""\
TDS_univariate divide=%s
  attributes='%s'
  filters='%s'
  options='%s'""" % (divi, att1, fil1, opt1))
  Canvas.cd(1); tds.draw(att1, fil1, opt1)
  Canvas.cd(2); tds.drawIter(att1, fil1, opt1)
  Canvas.cd(3); tds.drawCDF(att1, fil1, opt1)
  Canvas.cd(4); tds.drawBoxPlot(att1, fil1, opt1)
  _tdsNoGarbage.append(tds)
  return tds


###############################################################
def getStdReplaces():
  """standart useful replaces, @xxx@ as .in autotools"""
  replaces = [
    ("@WHOAMI@", getpass.getuser()),
    ("@DATE@", datetime.now().strftime("%Y/%m/%d")),
    ("@TIME@", datetime.now().strftime("%H:%M:%S")),
    ("@DATETIME@", datetime.now().strftime("%y%m%d_%H%M%S")),
    ("@BASEFILES@", "\n  ".join(getBaseFiles())),
  ]
  return replaces

###############################################################
def getBaseFiles():
  """list of useful mandatory files in uranie etude directory"""
  res = r"README.txt Doxyfile".split() #mandatories
  return res

###############################################################
def listOnlyDirs(rootDir):
  """return list of expanded directories names in directory"""
  res = os.listdir(rootDir)
  res = [os.path.join(rootDir, i) for i in res]
  return [i for i in res if os.path.isdir(i)]

###############################################################
def isPresentCMakeLists(aDir):
  """test if CMakeLists.txt exists in directory"""
  return os.path.isfile(os.path.join(aDir, "CMakeLists.txt"))

###############################################################
def getFixedUsefulDirs(rootDir):
  """unconditionaly fixed useful"""
  fixed = [os.path.join(rootDir, i) for i in "data html macros".split()]
  return fixed

###############################################################
def getOtherUsefulDirsForCpack(rootDir):
  """
  unconditionaly fixed useful are not in other useful directories
  get user subdirs with CMakeLists.txt as useful
  """
  fixed = getFixedUsefulDirs(rootDir)

  res = listOnlyDirs(rootDir)
  res = [i for i in res if i not in fixed] #to avoid test isPresentCMakeLists on fixed
  #print "OtherUsefulDirs:", res
  res = [i for i in res if isPresentCMakeLists(i)]
  return res

###############################################################
def getUsefulDirsForCpack(rootDir):
  #unconditionaly useful
  fixed = getFixedUsefulDirs(rootDir)
  
  res = listOnlyDirs(rootDir)
  res = [i for i in res if i not in fixed] #to avoid test isPresentCMakeLists on fixed
  #print "UsefulDirs:", res
  res = [i for i in res if isPresentCMakeLists(i)]
  return fixed + res #extend

###############################################################
def getUselessDirsForCpack(rootDir):
  """list of useless directories in uranie etude directory, not packaging saved"""
  #unconditionaly useless
  fixed = [os.path.join(rootDir, i) for i in "CMakeFiles _CPack_Packages .git".split()]
  
  res = listOnlyDirs(rootDir)
  res.extend(fixed) #useless even not existing yet

  useful = getUsefulDirsForCpack(rootDir)
  res = [i for i in res if i not in useful]
  return res

###############################################################
def _getAllFiles(rootDir, pattern="*"):
  """directories recursive"""
  matches = []
  for root, dirnames, filenames in os.walk(rootDir):
      if '.git' in root: continue
      for filename in fnmatch.filter(filenames, pattern):
          matches.append(os.path.join(root, filename))
  return matches

###############################################################
def getUselessFilesForCpack(rootDir, useful):
  """list of usefless files in uranie etude directory, not packaging saved"""
  res = """\
*~
CMakeCache.txt
CMakeLists.txt.user*
CPackConfig.cmake
CPackSourceConfig.cmake
install_manifest.txt
cmake_install.cmake
CTestTestfile.cmake
Makefile
""".split()

  res = [os.path.join(rootDir, i) for i in res] #unconditionaly useless
  res.extend(_getAllFiles(rootDir))
  #res = [i for i in res if os.path.isfile(i)]
  res = [i for i in res if i not in useful]
  return res

###############################################################
def getIgnoreFilesForCpack(rootDir, useful, Verbose=False):
  """list of usefless directories in uranie etude directory, not packaging saved"""
  dirs = getUselessDirsForCpack(rootDir)
  files = getUselessFilesForCpack(rootDir, useful)
  res = dirs + files
  keepFiles = "uranieGui.xml CMakeLists.txt".split() #want them in cpack
  res = [i for i in res if os.path.basename(i) not in keepFiles]
  #keep all files in subdirs other else data html macros with CMakeLists.txt
  otherDirsToKeep = getOtherUsefulDirsForCpack(rootDir)

  res = _filterOtherDirsToKeep(res, otherDirsToKeep)
  
  if Verbose: 
    print("UselessDirsForCpack\n%s" % PP.pformat(dirs))
    print("UselessFilesForCpack\n%s" % PP.pformat(files))
    print("usefulFilesForCpack\n%s" % PP.pformat(useful))
    print("IgnoreFilesForCpack\n%s" % PP.pformat(res))
  return res

def _filterOtherDirsToKeep(listFiles, otherDirsToKeep):
  """
  listFiles files will cpack ignored
  remove files to keep in .tgz from listFiles
  """
  res = []
  for afile in listFiles:
    tokeep = False #to keep in .tgz
    for dirkeep in otherDirsToKeep:
      if dirkeep in afile: 
        tokeep = True
        break
    if not tokeep: 
      res.append(afile)
    
  #print "filterOtherDirsToKeep\n",listFiles,"\n->\n",res
  return res

###############################################################
def execReplaces(aStr, replaces=[]):
  """append standart useful replaces to users replaces"""
  #users replaces executed before...
  replacesCurrent = [i for i in replaces]
  replacesCurrent.extend(getStdReplaces())
  res = str(aStr)
  for ini, fin in replacesCurrent:
    #print "ini, fin", ini, fin
    res = res.replace(ini, fin)
  return res

###############################################################
def getFilePatterns(name, replaces=[]):
  try:
    aStr = _filePatterns[str(name)]
  except:
    logger.error("no pattern for '%s', patterns are:\n%s" % (name, PP.pformat(sorted(_filePatterns.keys()))))
    return ""
  return execReplaces(aStr, replaces)

###############################################################
def getPatternKeys():
  return sorted(_filePatterns.keys())

###############################################################
# utilities for .dat file headers
###############################################################

###############################################################
def filterColumsNamesForUranie(header):
  """
  assume inexisting header array with empty strings ''
  lenght of #COLUMN_NAMES
  fill incomplete array with ''
  """
  def get_headerArray(name):
    try:
      nb = len(header["COLUMN_NAMES"])
    except:
      nb = 0
    res = ['' for i in range(nb)]
    try:
      ini = header[name]
      for i in range(nb): res[i] = ini[i]
      return res
    except:
      pass
    return res

  res = []
  #fmt = "%s{%i, '%s', '%s'}"
  fmt = "#COLUMN: %s | %i | %s | %s"
  try:
    cn = header["COLUMN_NAMES"]
    ct = get_headerArray("COLUMN_TITLES")
    cu = get_headerArray("COLUMN_UNITS")
    nb = len(cn)
    tags = [isUranieColumn(i) for i in cn]
    res = ['' for i in range(nb)]
    res = [fmt % (cn[i], i, ct[i], cu[i]) for i in range(nb) if tags[i]]
    res.extend([fmt % (cn[i], i, ct[i], cu[i]) for i in range(nb) if not tags[i]])
  except:
    pass
  return res

###############################################################
def isUranieColumn(name):
  uranieTags = 'n__iter__'.split() #may be some, later
  for tag in uranieTags:
    if tag in name: return True
  return False

###############################################################
def getHeaderContentsFileDat(filenameIni):
  """
  header description of .dat uranie/salome/paraview
  returns a dict with keys NAME,DATE...
  message if inexisting file
  """
  filename = os.path.expandvars(str(filenameIni))
  if os.path.isfile(filename):
    header = getHeaderFileDat(filename)
    return header
  else:
    logger.warning("inexisting file '%s'" % filename)
    return {}

###############################################################
def getHeaderFileDat(filename, Verbose=False):
  """
  header description of .dat uranie/salome/paraview
  returns a dict with keys #NAME, #DATE...
  """
  head = []
  with open(filename, "r") as f: 
    try:
      for i in xrange(20): #20 first lines
        head.append(f.next().strip())
    except: #StopIteration if shorter file
      pass
  #accept lines beginning '#', skip lines '#' alone
  head = [i for i in head if ("#" in i) and (i != "#")]
  res = {}
  for i in head:
    ok, tag, value = filterHeaderDat(i)
    if ok: 
      res[tag] = value
    else:
      if Verbose: 
        print("problem header file .dat:\n  '%s':\n  %s" % (filename, value))
      pass
  return res

###############################################################
def filterHeaderDat(line):
  """
  extract (key, value) from line from file .dat uranie/salome/paraview
  """
  """
  examples:
  "#COLUMN_UNITS: cm | cm |..." returns ('COLUMN_UNITS', ['cm','cm',...])
  "#NAME: iris" returns ('NAME', 'iris')
  """
  sp = line.split(": ")
  if sp[0][0] != "#":
    return (False, None, "line have to begin by '#': '%s'" % line)
  try:
    r1 = sp[0][1:]
    r2 = sp[1]
    if "COLUMN_" in r1:
      tr2 = [i.strip() for i in r2.split("|")]
      return (True, r1, tr2)
    else:
      return (True, r1, r2.strip())
  except:
    #logger.error("problem with '%s'" % line)
    return (False, None, "error in line: '%s'" % line)

 



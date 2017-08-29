"""
Warping of well logs from Teapot Dome survey.
Author: Loralee Wheeler, Colorado School of Mines
Author: Dave Hale, Colorado School of Mines
Version: 2015.01.30
"""

from java.awt import *
from java.io import *
from java.nio import *
from java.lang import *
from java.util import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mesh import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from warpt import *

##########################################################################

wlw = WellLogWarping()
wells = None
curves  = ["v", "p", "d"] # velocity, porosity, density
weights = [2.0, 1.0, 2.0] # set any weight to zero to exclude that log type
epow    = [0.25,0.25,0.125] # power of norm for alignment errors 

ms = 500 # maximum shift
fnull = -999.2500 # placeholder for missing data

##########################################################################

def main(args):
  global logs,sz
  sz,logs = getLogs(curves)
  goShifts()
  goWarping()
  #goErrors()

"""
Performs simultaneous correlation of many logs. 
Plots logs before and after alignment.
"""
def goShifts():
  wl = logs
  nz,nc,nl = len(wl[0][0]),len(wl[0]),len(wl)
  print(nz, nc, nl)
  wlw.setPowError(epow)
  wlw.setMaxShift(ms)
  s = wlw.findShifts(weights, wl)
  ##s = mul(1000*sz.delta,s) 
  s = mul(sz.delta,s) 
  fs = zerofloat(nz,nl)
  freplace = 2.0
  #fclips = (2.0,6.0)
  fclips = (40.0,190.0) # we use min-max of the DT values in the file vlogs.txt is just the name 
  #cblabel = "Velocity (km/s)" #velocity
  cblabel = "DT (.US/F)" 
  for i,c in enumerate(curves):
    if (weights[i]>0):
      for j in range(nl):
        fs[j] = wl[j][i]
      gs = wlw.applyShifts(fs,s)
      if c=="d":
        freplace = 1.0 # step for colourbar
        #fclips = (2.0,2.8)  
        fclips = (10.0, 200) #colourbar axis min-max
        #cblabel = "Density (g/cc)"
	cblabel = "Gamma Ray (GAPI)"
      if c=="p":
        freplace = 1.0 # step for colourbar 
        fclips = (5.0, 170) #colourbar axis min-max
        cblabel = "NPHI (.V/V)"
      """
      if c=="p":
        freplace = 0.0
        fclips = (1.0,100.0)
        cblabel = "Porosity (%)"
      fs = wlw.replaceNulls(fs,freplace)
      gs = wlw.replaceNulls(gs,freplace)
      if c=="p":
        fs = mul(fs,100)
        gs = mul(gs,100)
      """
      sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
      sp.setSize(500,900)  #size of plot
      #sp.setVLabel("Depth (km)")
      sp.setVLabel("Depth (m)")
      sp.setHLabel("Log index")
      sp.addColorBar(cblabel)
      sp.plotPanel.setColorBarWidthMinimum(90)
      pv = sp.addPixels(sz,Sampling(nl,1,1),fs)
      pv.setInterpolation(PixelsView.Interpolation.NEAREST)
      pv.setColorModel(ajet)
      pv.setClips(fclips[0],fclips[1])

      sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
      sp.setSize(500,900)
      sp.setVLabel("Relative geologic time")
      sp.setHLabel("Log index")
      sp.addColorBar(cblabel)
      sp.plotPanel.setColorBarWidthMinimum(90)
      pv = sp.addPixels(sz,Sampling(nl,1,1),gs)
      pv.setInterpolation(PixelsView.Interpolation.NEAREST)
      pv.setColorModel(ajet)
      pv.setClips(fclips[0],fclips[1])

"""
Performs pairwise correlation between specified log pairs.
Plots log pairs before and after alignment.
"""
def goWarping():
  # define pairs of logs to compare below
  #pairs = [(0,1),(2,4),(4,11)]
  pairs = [(0,1),(2,4),(4,6)] 
  wl = logs
  wlw.setPowError(epow)
  wlw.setMaxShift(ms)
  freplace = 2.0
  #fclips = (2.0,6.0)
  fclips = (40.0,190.0) # we use min-max of the DT values in the file vlogs.txt is just the name 
  #cblabel = "Velocity (km/s)" #velocity
  cblabel = "DT (.US/F)"
  for ic,c in enumerate(curves):
    if c=="d":
      freplace = 2.0
      #fclips = (1.8,3.4)
      fclips = (10,140)
      cblabel = "Gamma Ray (.GAPI)"
      ltype = "Gamma Ray"
    if c=="p":
      #freplace = 0.0
      #fclips = (0.0,0.6)
      #cblabel = "Porosity"
      #ltype = "Porosity "
      freplace = 0.0
      fclips = (5.0,140.0)
      cblabel = "NPHI"
      ltype = "NPNI "
    if (weights[ic]>0):
      for pair in pairs:
        jf,jg = pair[0],pair[1]
        fi,gj = wl[jf][ic],wl[jg][ic]
        if (wlw.wellNotNull(fi) and wlw.wellNotNull(gj)):
          e = wlw.computeErrors(ic,fi,gj)
          nl,nk = len(e[0]),len(e)
          d = wlw.accumulateErrors(e)
          kl = wlw.findWarping(d)
          fk,gk = wlw.applyWarping(kl,fi,gj)
          fi = wlw.replaceNulls(fi,freplace)
          gj = wlw.replaceNulls(gj,freplace)
          fk = wlw.replaceNulls(fk,freplace)
          gk = wlw.replaceNulls(gk,freplace)
          ii,ff = removeZeros(fi)
          jj,gg = removeZeros(gj)
          ki,kf = removeZeros(fk)
          kj,kg = removeZeros(gk)
          for ij in range(len(ff)):
            if ff[ij] == -2.0:
              ff[ij] = 3.0
          for ij in range(len(gg)):
            if gg[ij] == -2.0:
              gg[ij] = 3.0
          for ij in range(len(kf)):
            if kf[ij] == -2.0:
              kf[ij] = 3.0
          for ij in range(len(kg)):
            if kg[ij] == -2.0:
              kg[ij] = 3.0
          sff = Sampling(len(ff),sz.getDelta(),sz.getFirst()+ii*sz.getDelta())
          sgg = Sampling(len(gg),sz.getDelta(),sz.getFirst()+jj*sz.getDelta())
          skf = Sampling(len(kf),0.5*sz.getDelta(),sz.getFirst()+ki*0.5*sz.getDelta())
          skg = Sampling(len(kg),0.5*sz.getDelta(),sz.getFirst()+kj*0.5*sz.getDelta())
          title = ltype+"("+str(jf)+","+str(jg)+")"
          if True:
            sp = SimplePlot()
            sp.setSize(750,500)
            sp.setTitle(title+" before")
            #sp.setHLabel("Depth (km)")
	    sp.setHLabel("Depth (m)")
            sp.setVLabel(cblabel)
            sp.setVLimits(fclips[0],fclips[1])
            sp.setHLimits(0,1.9)
            pv = sp.addPoints(sff,ff)
            pv.setLineColor(Color.BLACK)
            pv = sp.addPoints(sgg,gg)
            pv.setLineColor(Color.RED)
          if True:
            sp = SimplePlot()
            sp.setSize(750,500)
            sp.setTitle(title+" after")
            #sp.setHLabel("Depth (km)")
            sp.setHLabel("Depth (m)")
            sp.setVLabel(cblabel)
            sp.setVLimits(fclips[0],fclips[1])
            sp.setHLimits(0,1.9)
            pv = sp.addPoints(skf,kf)
            pv.setLineColor(Color.BLACK)
            pv = sp.addPoints(skg,kg)
            pv.setLineColor(Color.RED)


"""
Computes alignment errors on the kl-coordinate system for specified log 
pairs and finds an optimal path.
Plots alignment error array with and without an optimal path.
"""
def goErrors():
  # define pairs of logs to compare below
  #pairs = [(0,1),(2,4),(4,11)] 
  pairs = [(0,1),(2,4),(4,6)] 
  wl = logs
  wlw.setPowError(epow)
  wlw.setMaxShift(ms)
  freplace = 2.0
  fclips = (2.0,6.0)
  cblabel = "Velocity (km/s)"
  ltype = "Velocity "
  for ic,c in enumerate(curves):
    if c=="d":
      freplace = 1.0
      fclips = (1.8,3.4)
      cblabel = "Density (g/cc)"
      ltype = "Density "
    if c=="p":
      freplace = 0.0
      fclips = (0.0,0.6)
      cblabel = "Porosity"
      ltype = "Porosity "
    if (weights[ic]>0):
      for pair in pairs:
        jf,jg = pair[0],pair[1]
        fi,gj = wl[jf][ic],wl[jg][ic]
        if (wlw.wellNotNull(fi) and wlw.wellNotNull(gj)):
          e = wlw.computeErrors(ic,fi,gj)
          wlw.interpolateOddErrors(e)
          nl,nk = len(e[0]),len(e)
          lmax = (nl-1)/2
          lmin = -lmax
          sl = Sampling(nl,1,lmin)
          sk = Sampling(nk,1,0)

          title = ltype+"("+str(jf)+","+str(jg)+")"
          # plots alignment errors 
          sp = SimplePlot()
          sp.setSize(750,500)
          sp.setTitle(title)
          sp.setHLabel("Depth index k")
          sp.setVLabel("Lag index l")
          sp.setHFormat("%5f")
          pv = sp.addPixels(sk,sl,transpose(e))
          pv.setInterpolation(PixelsView.Interpolation.NEAREST)
          d = wlw.accumulateErrors(e)
          wlw.interpolateOddErrors(d)

          # plots alignment errors with optimal path
          sp = SimplePlot()
          sp.setSize(750,500)
          sp.setTitle(title)
          sp.setHLabel("Depth index k")
          sp.setVLabel("Lag index l")
          sp.setHFormat("%5f")
          pv = sp.addPixels(sk,sl,transpose(e))
          pv.setInterpolation(PixelsView.Interpolation.NEAREST)
          d = wlw.accumulateErrors(e)
          wlw.interpolateOddErrors(d)
          kw,lw = wlw.findWarping(d)
          kw = wlw.toFloat(kw)
          lw = wlw.toFloat(lw)
          pv = sp.addPoints(kw,lw)
          pv.setLineColor(Color.WHITE)

          
          # plot accumulated alignment errors
          sp = SimplePlot()
          sp.setSize(750,500)
          sp.setTitle(title)
          sp.setHLabel("Depth index k")
          sp.setVLabel("Lag index l")
          sp.setHFormat("%5f")
          pv = sp.addPixels(sk,sl,transpose(d))
          #pv.setColorModel(cjet)
          pv.setInterpolation(PixelsView.Interpolation.NEAREST)
          pv.setPercentiles(0,90)
          kw,lw = wlw.findWarping(d)
          kw = wlw.toFloat(kw)
          lw = wlw.toFloat(lw)
          pv = sp.addPoints(kw,lw)
          pv.setLineColor(Color.WHITE)
          

#############################################################################
# utilities

def getLogs(curves):
  nc = len(curves)
  print("Number of curves:", nc)

  fileName = "data/"+curves[0]+"logs.txt"
  ifile = open(fileName,'r+')
  lines = ifile.readlines()
  nz = len(lines)
  ifile.close()
  ifile = open(fileName,'r+')  #open for read/write
  c = 0  #counts how many lines do not contain well log data
  start = False
  while (start!=True):
    line = ifile.readline()
    c += 1
    if line == "":
      print('End of file')
      break
    elif line[0] == '~':
      start = True

  # the file is still open and the next line contains data
  line  = ifile.readline()
  wdata = line.split('\t')
  nl = len(wdata)-1  #number of logs

  logs = zerofloat(nz-c,nc,nl-1) #initialising
  depth = zerodouble(nz-c)
  #maybe ifile.close()
  for ic,cv in enumerate(curves):
    fileName = "data/"+cv+"logs.txt"
    ifile = open(fileName,'r+')
    c = 0
    start = False
    while (start!=True):
      line = ifile.readline()
      c += 1
      if line == "":
        print('End of file')
        break
      elif line[0] == '~':  # we just want this as the first character to start
        start = True

    line = ifile.readline()
    wdata = line.split('\t')  # of course, you should use the same format :p

    temp_data = wdata[0].split(' ')
    len1 = len(temp_data)
    if len1 > 1:
        depth[0] = float(temp_data[len1-1])
    else:
        depth[0] = float(wdata[0])

    for l in range(1,nl):	
      logs[l-1][ic][0] = float(wdata[l])  #something like a 3D array
    i = 1
    start = False
    while (start!=True):
      line = ifile.readline()
      if line == "":
        start = True
      else:
        wdata = line.split('\t')

        temp_data = wdata[0].split(' ')
	len1 = len(temp_data)       
        if len1 > 1:
            depth[i] = float(temp_data[len1-1])
        else:
            depth[i] = float(wdata[0])

        for l in range(1,nl):
          logs[l-1][ic][i] = float(wdata[l])
        i += 1
    ifile.close()
  zs = Sampling(depth)  #mapping from sample value to function value
  sz = Sampling(zs.count,zs.delta,zs.first)
  ##sz = Sampling(zs.count,zs.delta*0.001,zs.first)	
  print("zs", zs)
  return sz,logs

def removeZeros(f):
  n = len(f)
  i = 0
  while f[i] == -2.0:
    i += 1
  fs = zerofloat(n)
  c = 0
  for j in range(i,n):
    fs[c] = f[j]
    if (f[j] == -2.0 and j+2 >= n) or (f[j] == -2.0 and f[j+1] == -2.0
        and f[j+2] == -2.0):
      break
    c += 1
  #ft = zerofloat(c)
  ft = copy(c,fs)
  return i,ft




#############################################################################
# graphics

cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)

#############################################################################
# Run the function main on the Swing thread
import sys
from javax.swing import *
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

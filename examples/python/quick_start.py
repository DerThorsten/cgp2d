import vigra
import cgp2d
import numpy
import opengm
from termcolor import colored
import sys
import gc

img 		= vigra.readImage('lena.bmp')
img 		= vigra.sampling.resize(img,shape=[100,100])

scale 		= 1.5
gradMag 	= vigra.filters.gaussianGradientMagnitude(img,scale)
labels,nSeg = vigra.analysis.watersheds(gradMag)
labels		= numpy.require(labels,dtype=numpy.uint64)



tgrid 	= cgp2d.TopologicalGrid(labels)
cgp  	= cgp2d.Cgp(tgrid)


img 		= vigra.sampling.resize(img,shape=cgp.shape)

#cgp2d.visualize(img,cgp)


cells1=cgp.cells(1)
cell1=cells1[0]

b = cell1.bounds
vb= cell1.boundView()

print colored('bounds     ','green') , numpy.array(b)
print colored('boundsView ','blue')  , vb

vb[0]=22

print colored('bounds     ','green') , numpy.array(b)
print colored('boundsView ','blue')  , vb


cgp = None
cells1=None
cell1=None
del cells1
del cell1
gc.collect()


#print colored('bounds     ','green') , numpy.array(b)
print colored('boundsView ','blue')  , vb




sys.exit(0)






for cellType in [1] :
	for cell in cgp.cells(cellType):

		print colored('cellType','green'),':',	colored(cellType,'red')

		print colored('label 	 ','green'),	cell.label
		print colored('bounds        ','green'),	numpy.array(cell.bounds)
		print colored('boundsView    ','blue'),	cell.boundView()



if False :
	for cellType in [0,1,2] :
		for cell in cgp.cells(cellType):

			print colored('cellType','green'),':',	colored(cellType,'red')

			print colored('label 	 ','green'),	cell.label
			print colored('bounds    ','green'),	numpy.array(cell.bounds)
			print colored('bounded by','green'),	numpy.array(cell.boundedBy)
			print colored('coordinats','green'),	numpy.array(cell.pointArray())




"""
numVar 					= cgp.numCells(2)
numSecondOrderFactors 	= cgp.numCells(1)
numHighOrderFactors		= cgp.numCells(0)

print numVar,numSecondOrderFactors,numHighOrderFactors



cells0 = cgp.cells(0)
cells1 = cgp.cells(1)
cells2 = cgp.cells(2)


for cell0 in cells0:
	print cell0.label
	cell1Bounds  = numpy.array(cell0.bounds)-1
	cell2Bounds  = set()
	for cell1Index cell1Bounds:
"""		
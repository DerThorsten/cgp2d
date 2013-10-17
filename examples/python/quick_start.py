import vigra
import cgp2d
import numpy
from termcolor import colored

img 		= vigra.readImage('lena.bmp')
img 		= vigra.sampling.resize(img,shape=[15,15])

scale 		= 2.5
gradMag 	= vigra.filters.gaussianGradientMagnitude(img,scale)
labels,nSeg = vigra.analysis.watersheds(gradMag)
labels		= numpy.require(labels,dtype=numpy.uint64)



tgrid 	= cgp2d.TopologicalGrid(labels)
cgp  	= cgp2d.Cgp(tgrid)


img 		= vigra.sampling.resize(img,shape=cgp.shape)

#cgp2d.visualize(img,cgp)




for cellType in [0,1,2] :
	for cell in cgp.cells(cellType):

		print colored('cellType','green'),':',	colored(cellType,'red')

		print colored('label 	 ','green'),	cell.label
		print colored('bounds    ','green'),	numpy.array(cell.bounds)
		print colored('bounded by','green'),	numpy.array(cell.boundedBy)
		print colored('coordinats','green'),	numpy.array(cell.pointArray())


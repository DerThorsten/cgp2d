import vigra
import cgp2d
import numpy
import opengm
from termcolor import colored
import sys
import gc
import pylab
import matplotlib.cm as cm
import scipy.ndimage

img = vigra.readImage('lena.bmp')#[0:100,0:200,0:3]

print img.shape


imgHist = cgp2d.jointColorHistogram(
	img,
	tuple([0.0  ]*3),
	tuple([255.0]*3),
	tuple([3]*3),
	1
)


sigma = [1.2]*2 + [1.2]*3

smmothed = scipy.ndimage.filters.gaussian_filter(imgHist,sigma=sigma,order=0,mode='constant',cval=0.0)



h0 = imgHist.reshape( [imgHist.shape[0],imgHist.shape[1],-1   ])
h1 = smmothed.reshape( [smmothed.shape[0],smmothed.shape[1],-1   ])


# show the histogram

for x in range(0,h0.shape[2],3):

	f = pylab.figure()

	f.add_subplot(1, 3, 1)  # this line outputs images side-by-side
	pylab.imshow(numpy.swapaxes(img/255.0,0,1),cmap=cm.Greys_r)

	hImg = h0[:,:,x:x+3].copy()
	hImg-=hImg.min()
	hImg/=hImg.max()
	f.add_subplot(1, 3, 2)  # this line outputs images side-by-side
	pylab.imshow(numpy.swapaxes(hImg,0,1),cmap=cm.Greys_r)


	hImg2 = h1[:,:,x:x+3].copy()
	hImg2-=hImg2.min()
	hImg2/=hImg2.max()
	f.add_subplot(1, 3, 3)  # this line outputs images side-by-side
	pylab.imshow(numpy.swapaxes(hImg2,0,1),cmap=cm.Greys_r)

	pylab.title('Double image')
	pylab.show()
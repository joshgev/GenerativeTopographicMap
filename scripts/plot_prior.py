import getopt
import matplotlib.pyplot as pyplot
import sys
import numpy
import math
import os

opts,args = getopt.getopt(sys.argv[1:],"st")

double = True
stacked = False


for o,v in opts:
	if o == "-s":
		double = False
	if o == "-t":
		stacked = True


def do_file(f):
	print "Processing file: "+f
	fin = open(f)

	type0=[]
	type1=[]
	for line in fin:
		elements = line.rstrip('\n').split()
		if double:
			if elements[0] == "0":
				type0.append([float(elements[1]),float(elements[2]),float(elements[3]) ])
			else:
				type1.append([float(elements[1]),float(elements[2]),float(elements[3]) ])
		else:
			type0.append([float(elements[1]),float(elements[2]),float(elements[3]) ])
	
	X0 = numpy.array([i[1] for i in type0])
	Y0 = numpy.array([i[2] for i in type0])
	Z0 = [i[0] for i in type0]
	max0 = max(Z0)
	max1 = 0
	if double:
	        X1 = numpy.array([i[1] for i in type1])
	        Y1 = numpy.array([i[2] for i in type1])
	        Z1 = [i[0] for i in type1]
		max1 = max(Z1)

	max_o = max(max0,max1)
	
	max_o *= 1.1

	X0 = X0.reshape((math.sqrt(X0.size),math.sqrt(X0.size)))
	Y0 = Y0.reshape((math.sqrt(Y0.size),math.sqrt(Y0.size)))	
	Z0 = numpy.array(Z0).reshape(len(X0),len(X0))

	if double:
		X1 = X1.reshape((math.sqrt(X1.size),math.sqrt(X1.size)))
	        Y1 = Y1.reshape((math.sqrt(Y1.size),math.sqrt(Y1.size)))
	        Z1 = numpy.array(Z1).reshape(len(X1),len(X1))
	
	if double:
		fig = pyplot.figure()
		if not stacked:
			sp1 = fig.add_subplot(211)
			contour0 = sp1.contourf(X0,Y0,Z0,cmap=pyplot.cm.Reds,levels=numpy.arange(0,max_o,max_o/10))
			cb0 = pyplot.colorbar(contour0,format='%0.4f')
	
			sp2 = fig.add_subplot(212)
			contour1 = sp2.contourf(X1,Y1,Z1,cmap=pyplot.cm.Blues,levels=numpy.arange(0,max_o,max_o/10))
			cb1 = pyplot.colorbar(contour1, format='%0.4f')
			fig.savefig("contours/"+f+".png")
		else:
			sp1 = fig.add_subplot(111)
                        contour0 = sp1.contourf(X0,Y0,Z0,cmap=pyplot.cm.Reds,levels=numpy.arange(0,max_o,max_o/10))
                        cb0 = pyplot.colorbar(contour0, format='%0.4f')

                        contour1 = sp1.contourf(X1,Y1,Z1,cmap=pyplot.cm.Blues,levels=numpy.arange(0,max_o,max_o/10),alpha = .5)
                        cb1 = pyplot.colorbar(contour1, format='%0.4f')
                        fig.savefig("contours/"+f+".png")


		pyplot.close()
	else:
		fig = pyplot.figure()
		sp = fig.add_subplot(111)
		contour0 = sp.contourf(X0,Y0,Z0,cmap=pyplot.cm.Reds,levels=numpy.arange(0,max_o,max_o/10))
                cb0 = pyplot.colorbar(contour0, format='%0.4f')

		fig.savefig("contours/"+f+".png")
		pyplot.close()

present = []

for path,dirnames,filenames in os.walk('contours/'):

	present += (filenames)

files = []
for f in args:

	if f+".png" in present:
		continue
	files.append(f)

for f in files:

	do_file(f)


#This script generates the data for a simple test of the GTM algorithm.
#A two gaussian processes live in a two dimensional space with centers at (-1,0) and (1,0).
#The standard deviation of the left and right distributions is specified on the command line.
#After generating N points from each distribution (specified on command line), a third dimension 
#is added to the data by draing from a uniform distribution over the interval [-1,1].
#Finally, the entire data set is rotated.


import numpy
import numpy.random as rand
import sys
from math import *

def rotation_matrix(phi, theta, psi):
	return numpy.array( 
	((cos(theta)*cos(psi), -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi), sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi)),
 	(cos(theta)*sin(psi),  cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi), sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi)), 
	(-sin(theta),      sin(phi)*cos(theta),                               cos(phi)*cos(theta) ) ) )

def gaussian2d(center, sigma):
	r = rand.normal(0,sigma)
	theta = rand.uniform() * 2*3.14159
	
	x = r * cos(theta)
	y = r * sin(theta)

	return (x+center[0],y+center[1])

def dist3(point):
	return sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2])


#usage: python script sigma_a sigma_b N

sigma_a = float(sys.argv[1])
sigma_b = float(sys.argv[2])
N = int(sys.argv[3])

theta = rand.uniform() * 2 * 3.14159
phi   = rand.uniform() * 2 * 3.14159
psi   = rand.uniform() * 2 * 3.14159

#rotation = rotation_matrix(theta,phi,psi)

fout = open ("gtm_latent.txt","w")

for i in range(N):

	distribution = "a" if rand.uniform() < .5 else "b"
	center = (-3,0) if distribution == "a" else (3,0)
	sigma  = sigma_a if distribution == "a" else sigma_b

	point = gaussian2d(center,sigma)

	z = rand.uniform(-1,1)
	
	fout.write(str(point[0]) + " "+str(point[1]) + "\n")
	point3d = (point[0],point[1],z)

	r = dist3(point3d)
	rotation = rotation_matrix(theta * r / .05,phi * r / .05 ,psi * r / .05)

	vector = numpy.array(point3d)
	final = numpy.dot(rotation  ,  vector)
	label = 1 if distribution == "a" else 0
	print str(label) +" "+str(final[0]) +" " +str(final[1])+" "+str(final[2])
		

	
	

fout.close()	

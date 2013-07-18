import matplotlib.pyplot as pyplot
import sys
import getopt

opts,args = getopt.getopt(sys.argv[1:],"x:y:sf")
x = 1
y = 2
multicolor = True
fixed_scale = False
for o,v in opts:
	if o == '-x':
		x = int(v)
	if o == '-y':
		y = int(v)
	if o== '-s':
		multicolor=False
	if o== '-f':
		fixed_scale=True
	
		
print "Using colums "+str(x)+" and "+str(y)
if multicolor:
	print "Using multicolor"

files = args
#print "Reading file: "+file
def do_file(f):
#for f in files:
	fin = open(f)
	x0 = []
	y0 = []
	x1 = []
	y1 = []
	for line in fin:
		line = line.rstrip("\n")
		elements = line.split()
		if multicolor:
			if (elements[0] == "0"):
				x0.append( float( elements[x] ) )
	
				y0.append( float( elements[y] ) )
			elif (elements[0] == "1"):
				x1.append( float( elements[x] ) )
				y1.append( float( elements[y] ) )
		else:
			x0.append( float( elements[x] ) )
			y0.append( float( elements[y] ) )
	
	fin.close()
	
	
	fig = pyplot.figure()
	ax1 = fig.add_subplot(111)
	if fixed_scale:
		ax1.axis([-1,1,-1,1])
	ax1.scatter(x0,y0,color='r')
	if multicolor:
		ax1.scatter(x1,y1,color='b',alpha=.5)
	if len(files) == 1:
		pyplot.show()
	else:	
		pyplot.savefig("plots/"+f+".png")
	fig.clf()
	pyplot.clf()
	pyplot.close()	

for f in files:
	print "plotting file: "+f
	do_file(f)

import random

out = open("pca_data.txt","w")
for i in range(0,1000):
	
	x = random.random();
	y = x + (random.random() - .5) / 10
	out.write(str(x) +" "+str(y)+"\n")
out.close()
	

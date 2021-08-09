import numpy as np

xx=np.zeros(3)
xx[0]=1


CC = np.zeros((3,3 ))
CC[1][1]=1
CC[2][2]=2
CC[0][0]=3

AA = np.zeros((3,3 ))
AA[1][1]=1
AA[2][2]=2
AA[0][0]=3

hilfe=sum(sum(AA[n][t]*CC[n][t]  for t in {0,1,2}) + xx[n] for n in {0,1,2} )

print(hilfe)
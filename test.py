import numpy as np
import matplotlib.pyplot as plt
import math
print("hack attempt")

data = np.genfromtxt('please', delimiter = '  ') #poner dos espacios para in.example9

# script para probar anexar vectores a matrices y trasponer 

matrix1 = np.array([[1,2,3],[4,5,6]])

#vector1 = matrix1[:,0:1]
#print(vector1)
'''
a=[1,2]
print(a)
b=np.matrix(a).T
#a=np.asarray(a)
#b= a.reshape(-1,1)
print(b)
c=np.concatenate((matrix1,b),1)
print(c)
'''

top=[0,0,0]
top=np.concatenate(([0],top))
#top=np.asarray(top)
#print(top)
#c = np.vstack((matrix1,top))


print(top)


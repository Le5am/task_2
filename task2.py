import numpy as np
import scipy.constants as sc
import scipy.special as ss
import wget
import matplotlib.pyplot as plt
import math as mh
import os
import xml.etree.ElementTree as ET

var = 2
dx = 1e+6  
url = 'https://jenyay.net/uploads/Student/Modelling/task_02.xml'
input_data = wget.download(url)
print(input_data)


tree = ET.parse('task_02.xml')
root = tree.getroot()

for elem in root.iter('variant'):
    if elem.get('number') == '2':
        D = float(elem.get('D'))
        fmin = float(elem.get('fmin'))
        fmax = float(elem.get('fmax'))
        
print (D, fmin, fmax)

f = np.arange(fmin, fmax, dx)
Lambda = sc.c/f
r = D/2
k = mh.pi*2/Lambda

#Bessele's functions
def hn(n,x): return ss.spherical_jn(n,x) + 1j*ss.spherical_yn(n,x)
def bn(n,x): return (x*ss.spherical_jn(n-1,x) - n*ss.spherical_jn(n,x))/(x*hn(n-1,x) - n*hn(n,x))
def an(n,x): return ss.spherical_jn(n,x)/hn(n,x)

arr_sum = [((-1)**n)*(n+0.5)*(bn(n,k*r) - an(n,k*r)) for n in range(1,20)]
sum = np.sum(arr_sum, axis=0)
sigma = (Lambda**2)/mh.pi * abs(sum)**2

if not os.path.isdir("results"):
    os.mkdir("results")
    
with open("./results/data.txt", "w") as file:
    for i in range(0, len(f)):
        file.write(f'{f[i]}    {sigma[i]}\n')

plt.plot(f/10e6, sigma)
plt.xlabel("f, МГц")
plt.ylabel("sigma, м^2")
plt.title('Задача 2')
plt.grid()
plt.show()
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 13:10:34 2018

@author: andy
"""

import numpy as np

N = 9981
alpha = 1e-4
pi = 0.8
T = 2

shift0 = np.zeros(N+1)
shift1 = np.zeros(N+1)

for i in range(N+1):
	if i<=99:
		shift0[i] = 1
		shift1[i] = 1
	elif i <= 900:
		shift0[i] = 1-alpha*np.log(i)**pi
		shift1[i] = 1-alpha*np.log(i)**pi
	else:
		shift0[i] = 1-alpha*np.log(i)**(2*pi)
		shift1[i] = 1-alpha*np.log(i)**(3*pi)
	

optimal0 = np.zeros((N+1,T))
optimal1 = np.zeros((N+1,T))

less1K_best = np.ones(10)*99
for k in range(1, 10):
	candidate = max(shift0[100*k]*100*k,less1K_best[k-1])
	if candidate > less1K_best[k]:
		less1K_best[k]=candidate
		
greater1K_best0 = np.ones(int(np.floor(N/1000.0))+1) * less1K_best[9]
greater1K_best1 = np.ones(int(np.floor(N/1000.0))+1) * less1K_best[9]

for k in range(1,len(greater1K_best0)):
	candidate = max(shift0[1000*k]*1000*k, greater1K_best0[k-1])
	if candidate > greater1K_best0[k]:
		greater1K_best0[k] = candidate
for k in range(1,len(greater1K_best1)):
	candidate = max(shift1[1000*k]*1000*k, greater1K_best1[k-1])
	if candidate > greater1K_best1[k]:
		greater1K_best1[k] = candidate		

for i in range(min(100,N+1)):
	optimal0[i,T-1] = i
	optimal1[i,T-1] = i

for i in range(1, min(10,int(np.floor(N/100))+1)):
	optimal0[i*100:min(i*100+99,N)+1,T-1] = less1K_best[i]
	optimal1[i*100:min(i*100+99,N)+1,T-1] = less1K_best[i]

for i in range(1, int(np.floor(N/1000))+1):
	optimal0[i*1000:min(i*1000+999,N)+1,T-1] = greater1K_best0[i]
	optimal1[i*1000:min(i*1000+999,N)+1,T-1] = greater1K_best1[i]
		
for t in np.arange(T-2, -1, -1):
	for i in range(N+1):
		bestone0 = 0
		bestone1 = 0
		for k in range(0, min(99,i)+1):
			candidate = shift0[k]*(k+optimal0[i-k,t+1])
			if candidate > bestone0:
				bestone0 = candidate
			if candidate > bestone1:
				bestone1 = candidate
		for k in range(1, min(9,int(np.floor(i/100.0)))+1):
			candidate = shift0[100*k]*(100*k+optimal1[i-100*k,t+1])
			if candidate > bestone0:
				bestone0 = candidate
			if candidate > bestone1:
				bestone1 = bestone1
		for k in range(1, int(np.floor(i/1000.0))+1):
			candidate = shift0[1000*k]*(1000*k+optimal1[i-1000*k,t+1])
			if candidate > bestone0:
				bestone0 = candidate
			candidate = shift1[1000*k]*(1000*k+optimal1[i-1000*k,t+1])
			if candidate > bestone1:
				bestone1 = candidate
		optimal0[i,t] = bestone0
		optimal1[i,t] = bestone1
		
print(optimal0[N,0])
		
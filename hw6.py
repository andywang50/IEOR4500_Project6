# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 13:10:34 2018

@author: andy
"""

import numpy as np

N = 9981
alpha = 1e-4
pi = 0.8
T = 10

shift0 = np.zeros(N+1)
shift1 = np.zeros(N+1)

execution0 = np.zeros((N+1,T))
execution1 = np.zeros((N+1,T))

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
t
for i in range(N+1):
	if i<=99:
		optimal0[i,T-1] = i
		optimal1[i,T-1] = i
		execution0[i,T-1] = i
		execution1[i,T-1] = i

	elif i <= 900:
		bestone = 99
		bestk = 99
		for k in range(1, min(9,int(np.floor(i/100.0)))+1):
			candidate = shift0[100*k]*100*k
			if candidate > bestone:
				bestone=candidate
				bestk=100*k
		optimal0[i,T-1] = bestone
		optimal1[i,T-1] = bestone
		execution0[i,T-1] = bestk
		execution1[i,T-1] = bestk

	else:
		bestone0 = 99
		bestone1 = 99
		bestk0 = 99
		bestk1 = 99
		for k in range(1, 10):
			candidate = shift0[100*k]*100*k
			if candidate > bestone0:
				bestone0=candidate
				bestk0 = 100*k
			if candidate > bestone1:
				bestone1 = candidate
				bestk1 = 100*k
		for k in range(1, int(np.floor(i/1000.0))+1):
			candidate = shift0[1000*k]*1000*k
			if candidate > bestone0:
				bestone0 = candidate
				bestk0 = 1000*k
			candidate = shift1[1000*k]*1000*k
			if candidate > bestone1:
				bestone1 = candidate
				bestk1=1000*k
		optimal0[i,T-1] = bestone0
		optimal1[i,T-1] = bestone1
		execution0[i,T-1] = bestk0
		execution1[i,T-1] = bestk1	
		
for t in np.arange(T-2, -1, -1):
	for i in range(N+1):
		bestone0 = 0
		bestone1 = 0
		bestk0 = 0
		bestk1 = 0
		for k in range(0, min(99,i)+1):
			candidate = shift0[k]*(k+optimal0[i-k,t+1])
			if candidate > bestone0:
				bestone0 = candidate
				bestk0=k
			if candidate > bestone1:
				bestone1 = candidate
				bestk1=k
		for k in range(1, min(9,int(np.floor(i/100.0)))+1):
			candidate = shift0[100*k]*(100*k+optimal1[i-100*k,t+1])
			if candidate > bestone0:
				bestone0 = candidate
				bestk0=100*k
			if candidate > bestone1:
				bestone1 = candidate
				bestk1=100*k
		for k in range(1, int(np.floor(i/1000.0))+1):
			candidate = shift0[1000*k]*(1000*k+optimal1[i-1000*k,t+1])
			if candidate > bestone0:
				bestone0 = candidate
				bestk0=1000*k
			candidate = shift1[1000*k]*(1000*k+optimal1[i-1000*k,t+1])
			if candidate > bestone1:
				bestone1 = candidate
				bestk1=1000*k
		optimal0[i,t] = bestone0
		optimal1[i,t] = bestone1
		execution0[i,t] = bestk0
		execution1[i,t] = bestk1

print("optimal value:", optimal0[N,0])
best_executions=np.zeros(T)
N_remaining = N
P = 1
prices = np.zeros(T)
large = 0
for t in range(0,T):
	if large == 0:
		bestk = int(execution0[N_remaining, t])
		P *= shift0[bestk]
	else:
		bestk = int(execution1[N_remaining, t])
		P *= shift1[bestk]

	if bestk >= 100:
		large = 1
	else:
		large = 0
	N_remaining -= bestk
	best_executions[t] = bestk
	prices[t] = P
print("executions:", best_executions)
print("prices:", prices)
print(np.dot(prices,best_executions))
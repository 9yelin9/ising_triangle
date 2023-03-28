# it.py : Ising model with triangle lattice

import os
num_thread = 1
os.environ['OMP_NUM_THREADS'] = str(num_thread)
os.environ['OPENBLAS_NUM_THREADS'] = str(num_thread)

import re
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-v', '--vector', type=float, default=None, help='<radious>\n')
parser.add_argument('-s', '--show', nargs='+', default=None, help='l                 : ShowObs\n'\
																  'o <fname> <dtype> : ShowObs\n'\
																  'e <fname>         : ShowErr\n')
args = parser.parse_args()                                                                     

class IsingTriangle:
	def __init__(self):
		self.DIM = 2
		self.NN  = 6

		self.obs_dict = {
			'e': {
				'col': 1,
				'full': 'Energy'
			},
			'm': {
				'col': 2,
				'full': 'Absolute magnetization'
			}
		}
		
	def CalcVector_(self, r, coef):
		return [n0 * r[0] + n1 * r[1] for n0 in coef[0] for n1 in coef[1]]

	def GenVector(self, rad=1):
		r_lat = np.array([[1, 0], [1/2, np.sqrt(3)/2]])
		r_sup = np.array([[2, np.sqrt(3)], [5/2, -np.sqrt(3)/2]])
		r_sub = []

		coef_lat = np.repeat([range(-4, 4)], repeats=self.DIM, axis=0)

		for r in self.CalcVector_(r_lat, coef_lat):
			if np.linalg.norm(r) < rad + 1e-6: r_sub.append(r)
		r_sub = np.array(r_sub)

		with open('input/vector.txt', 'w') as f:
			f.write('%20.16f\n\n' % rad)
			for r in [r_lat, r_sup, r_sub]:
				for ri in r: f.write('%s\n' % ''.join(['%20.16f' % rii for rii in ri]))
				f.write('\n')

		return r_lat, r_sup, r_sub
				
	def ShowLat(self, rad=1):
		r_lat, r_sup, r_sub = self.GenVector(rad=rad)

		coef_sup = np.repeat([range(-1, 2)], repeats=self.DIM, axis=0)

		fig, ax = plt.subplots(constrained_layout=True)

		ax.arrow(0, 0, *r_sup[0])
		ax.arrow(0, 0, *r_sup[1])

		for r in self.CalcVector_(r_sup, coef_sup):
			ax.scatter(r[0] + r_sub[:, 0], r[1] + r_sub[:, 1])

		plt.show()

	def DrawObs(self, fname, dtype, data, ax):
		ax.plot(data[:, 0], data[:, self.obs_dict[dtype]['col']], marker='.', label=fname)
		ax.set_title(self.obs_dict[dtype]['full'])
		ax.grid()
		ax.legend()

	def ShowObs(self, fname, dtype):
		fname = fname.split(',')
		dtype = dtype.split(',')

		data = []
		for fn in fname:
			with open('output/%s.txt' % fn, 'r') as f: data.append(np.genfromtxt(f))

		fig, ax = plt.subplots(1, len(dtype), figsize=(2+3*len(dtype), 4), constrained_layout=True)

		if len(dtype) > 1:
			for d, a in zip(dtype, ax):
				for fn, da in zip(fname, data): self.DrawObs(fn, d, da, a)
		else:
			for fn, da in zip(fname, data): self.DrawObs(fn, dtype[0], da, ax)

		fig.supxlabel('Temperature')
		fig.savefig('diagram/obs_%s_%s.png' % (''.join(fname), ''.join(dtype)))
		plt.show()

	def ShowErr(self, fname):
		fname = fname.split(',')

		fig, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)

		with open('output/test.txt', 'r') as f: data0 = np.genfromtxt(f)

		for fn in fname:
			with open('output/%s.txt' % fn, 'r') as f: data = np.genfromtxt(f)
			ax.plot(data0[:, 0], np.abs(data0[:, 1] - data[:, 1]), marker='.', label=fn)

		ax.grid()
		ax.legend()

		fig.supxlabel('Temperature')
		fig.supylabel(r'$|E_\mathrm{test} - E|$')
		fig.suptitle('Energy error')
		fig.savefig('diagram/err_%s.png' % ''.join(fname))
		plt.show()
###

it = IsingTriangle()
if args.vector: it.GenVector(rad=args.vector)
if args.show:
	if   args.show[0] == 'l': it.ShowLat()
	elif args.show[0] == 'o': it.ShowObs(*args.show[1:])
	elif args.show[0] == 'e': it.ShowErr(*args.show[1:])

#!~/.pyenv/shims/python3
# -*- coding: utf-8 -*-
#****************************************************************************#
# @Author: Ming Yang
# @Email:  ykevin09@gmail.com
# @Date:   2020-12-04 19:21:52
# @Last Modified by:   Ming Yang
# @Last Modified time: 2021-06-29 10:57:48
# @Function:
#****************************************************************************#
from math import gamma
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

plt.style.use(['classic', 'physics'])


def beta(x, a, b):
	B = gamma(a)*gamma(b)/gamma(a+b)
	res = integrate.quad(lambda t: (t**(a-1.))*((1.-t)**(b-1.)), 0., x)
	return 1/B*res[0]


def cvd_to_psd():
	""" CVD to PSD """
	infile = 'Ottawa-F65/Ottawa-F65_CMD.dat'
	d_e = np.loadtxt(infile, skiprows=1, delimiter=',')
	dmin, dmax = 0.00011, 0.00011*4
	dr = lambda x: (x - dmin)/(dmax - dmin)
	a, b = 2, 4 # 2, 4
	xx = np.linspace(0, 1, 100)
	yy = []
	for x in xx:
		yy.append(beta(x, a, b)*100)
	# PSD
	infile = 'spl/ClassGradingCurveFile.out'
	infile = 'spl/ClassSizeDistributionFile.out'
	d_p = np.loadtxt(infile, skiprows=1)
	d_V = []
	tmp = 0.
	for i, dp in enumerate(d_p):
		if i == 0:
			tmp = 1./6*(dp[0]**3)*dp[1]
			d_V.append(tmp)
		else:
			tmp = tmp + 1./6*(dp[0]**3)*(dp[1])
			d_V.append(tmp)
	#---- Make the plot ----
	plt.figure(0)
	plt.plot(dr(d_e[:, 0]*1e-3), d_e[:, 1], 'b-o', fillstyle='none', label=r'Ottawa-F65')
	plt.plot(xx, yy, 'g-.')
	plt.plot(dr(d_p[:, 0]), d_V/d_V[-1]*100, 'r--D', label=r'DEM')
	plt.xlabel(r'$d_{\text{r}}~(-)$')
	plt.ylabel(r'CVD~(\%)')
	plt.xlim(0, 1)
	plt.xticks(np.arange(0, 1.01, 0.2))
	plt.ylim(0, 100)
	plt.yticks(np.arange(0, 100.1, 20))
	plt.legend(loc='lower right', bbox_to_anchor=[1.015, -0.02])
	outfile = 'figs/DEM_Ottawa-F65_CVD_dr.eps'
	plt.savefig(outfile)
	#
	plt.figure(1)
	plt.plot(d_e[:, 0], d_e[:, 1], 'b-o', fillstyle='none', label=r'Ottawa-F65')
	plt.plot(d_p[:, 0]*1000, d_V/d_V[-1]*100, 'r--D', label=r'DEM')
	plt.xscale('log')
	plt.xlabel(r'$D~(\si{mm})$')
	plt.ylabel(r'CVD~(\%)')
	plt.xlim(0.01, 10)
	plt.ylim(0, 105)
	plt.legend(loc='lower right', bbox_to_anchor=[1.015, -0.02])
	outfile = 'figs/DEM_Ottawa-F65_CVD_D.eps'
	plt.savefig(outfile)
	#
	plt.show()


def spl_to_data():
	""" Sample info to data """
	# Number of atoms 
	infile = 'spl/ClassSizeDistributionFile.out'
	d_ = np.loadtxt(infile, skiprows=1)
	n_atom = int(np.sum(d_[:, 1]))
	n_atom_type = 2
	atom_type = 1
	rho = 2650.
	# SampleBoxSize
	infile = 'spl/SampleBoxSize.spl'
	d_ = np.loadtxt(infile, skiprows=1)
	xlo, ylo, zlo = d_[0, 0], d_[0, 1], d_[0, 2]
	xhi, yhi, zhi = d_[1, 0], d_[1, 1], d_[1, 2]
	#
	outfile = 'Ottawa-F65.spl'
	f_out = open(outfile, 'w')
	f_out.write('LIGGGHTS data file\n')
	f_out.write('\n')
	f_out.write(f' {n_atom} atoms\n')
	f_out.write('\n')
	f_out.write(f' {n_atom_type} atom types\n')
	f_out.write('\n')
	f_out.write(f' {xlo} {xhi} xlo xhi\n')
	f_out.write(f' {ylo} {yhi} ylo yhi\n')
	f_out.write(f' {zlo} {zhi} zlo zhi\n')
	f_out.write('\n')
	f_out.write('Atoms\n')
	f_out.write('\n')
	#
	infile = 'spl/sample.spl'
	f_in = open(infile, 'r')
	next(f_in)
	i_atom = 1
	for line in f_in:
		l_ = line.split(' ')
		if l_[0] == 'sphere':
			D, x, y, z = float(l_[2])*2, float(l_[3]), float(l_[4]), float(l_[5])
			f_out.write(f'{i_atom} {atom_type} {D} {rho} {x} {y} {z}\n')
			i_atom = i_atom + 1
	#
	f_out.close()
	f_in.close()


def main():
	""" Main function """
	# cvd_to_psd()
	spl_to_data()

if __name__ == '__main__':
	main()
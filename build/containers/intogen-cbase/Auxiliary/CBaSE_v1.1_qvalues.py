#!/usr/bin/env python

# ************************************************************
# * Cancer Bayesian Selection Estimation (CBaSE):			 *
# * Code accompanying Weghorn & Sunyaev, Nat. Genet. (2017). *
# *															 *
# * Author:		Donate Weghorn								 *
# *															 *
# * Copyright:	(C) 2017-2019 Donate Weghorn				 *
# *															 *
# * License:		Public Domain							 *
# *															 *
# * Version:	1.1											 *
# ************************************************************


# This version incorporates the following changes:
# - for sake of computation, this version skips the negative selection analysis;
# - this version changes the output format accordingly to display positive selection analysis;
# - add a selection clause for statistical model of choice based on the total mutation burden;
# - add an implementation of the Benjamini-Hochberg FDR method to compute q-values.


from scipy.optimize import minimize
from statsmodels.stats.multitest import multipletests
import scipy.special as sp
import sys
import math
import sys
import numpy as np
from statsmodels.stats.multitest import multipletests
import scipy.stats as st
import itertools as it
import subprocess
import scipy.special as sp
import mpmath as mp
import cgi
import cgitb
import argparse
import glob

# ************************************ FUNCTION DEFINITIONS *************************************


def compute_p_values(p, genes, aux):
	[modC, simC, runC] = aux
	if simC == 0:
		runC = 1

	if [1, 2].count(modC):
		a, b = p
	elif [3, 4].count(modC):
		a, b, t, w = p
	elif [5, 6].count(modC):
		a, b, g, d, w = p

	if modC == 1:
		# *************** lambda ~ Gamma:
		def pofs(s, L):
			return np.exp(s * np.log(L * b) + (-s - a) * np.log(1 + L * b) + sp.gammaln(s + a) - sp.gammaln(s + 1) - sp.gammaln(a))

		def pofx_given_s(x, s, L, r, thr):
			return np.exp(
				x * np.log(r) + (s + x) * np.log(L * b) + (-s - x - a) * np.log(1 + L * (1 + r) * b) + sp.gammaln(
					s + x + a) - sp.gammaln(s + 1) - sp.gammaln(x + 1) - sp.gammaln(a)) / pofs(s, L)
	elif modC == 2:
		# *************** lambda ~ IG:
		def pofs(s, L, thr):
			if thr:
				return 2. * mp.exp(
					((s + a) / 2.) * math.log(L * b) + mp.log(mp.besselk(-s + a, 2 * np.sqrt(L * b))) - sp.gammaln(
						s + 1) - sp.gammaln(a))
			else:
				return 2. * math.exp(
					((s + a) / 2.) * math.log(L * b) + np.log(sp.kv(-s + a, 2 * math.sqrt(L * b))) - sp.gammaln(
						s + 1) - sp.gammaln(a))

		def pofx_given_s(x, s, L, r, thr):
			if thr:
				return mp.exp(np.log(2) + (s + x) * np.log(L) + x * np.log(r) + (1 / 2. * (-s - x + a)) * np.log(
					(L * (1 + r)) / b) + a * np.log(b) + mp.log(
					mp.besselk(s + x - a, 2 * math.sqrt(L * (1 + r) * b))) - sp.gammaln(s + 1) - sp.gammaln(
					x + 1) - sp.gammaln(a)) / pofs(s, L, thr)
			else:
				return np.exp(np.log(2) + (s + x) * np.log(L) + x * np.log(r) + (1 / 2. * (-s - x + a)) * np.log(
					(L * (1 + r)) / b) + a * np.log(b) + np.log(
					sp.kv(s + x - a, 2 * math.sqrt(L * (1 + r) * b))) - sp.gammaln(s + 1) - sp.gammaln(
					x + 1) - sp.gammaln(a)) / pofs(s, L, thr)
	elif modC == 3:
		# *************** lambda ~ w * Exp + (1-w) * Gamma:
		def pofs(s, L):
			return np.exp(np.log(w) + s * np.log(L) + np.log(t) + (-1 - s) * np.log(L + t)) + np.exp(
				np.log(1. - w) + s * np.log(L * b) + (-s - a) * np.log(1 + L * b) + sp.gammaln(s + a) - sp.gammaln(
					s + 1) - sp.gammaln(a))

		def pofx_given_s(x, s, L, r, thr):
			return (np.exp(np.log(w) + (s + x) * np.log(L) + x * np.log(r) + np.log(t) + (-1 - s - x) * np.log(
				L + L * r + t) + sp.gammaln(1 + s + x) - sp.gammaln(s + 1) - sp.gammaln(x + 1)) + np.exp(
				np.log(1 - w) + x * np.log(r) + (s + x) * np.log(L * b) + (-s - x - a) * np.log(
					1 + L * (1 + r) * b) + sp.gammaln(s + x + a) - sp.gammaln(s + 1) - sp.gammaln(x + 1) - sp.gammaln(
					a))) / pofs(s, L)
	elif modC == 4:
		# *************** lambda ~ w * Exp + (1-w) * InvGamma:
		def pofs(s, L, thr):
			if thr:
				return (w * t * mp.exp(s * np.log(L) + (-1 - s) * np.log(L + t))) + mp.exp(
					np.log(1. - w) + np.log(2.) + ((s + a) / 2.) * np.log(L * b) + mp.log(
						mp.besselk(-s + a, 2 * math.sqrt(L * b))) - sp.gammaln(s + 1) - sp.gammaln(a))
			else:
				return (w * L ** s * t * (L + t) ** (-1 - s)) + np.exp(
					np.log(1. - w) + np.log(2.) + ((s + a) / 2.) * np.log(L * b) + np.log(
						sp.kv(-s + a, 2 * math.sqrt(L * b))) - sp.gammaln(s + 1) - sp.gammaln(a))

		def pofx_given_s(x, s, L, r, thr):
			if thr:
				return (np.exp(np.log(w) + (s + x) * np.log(L) + x * np.log(r) + np.log(t) + (-1 - s - x) * np.log(
					L + L * r + t) + sp.gammaln(1 + s + x) - sp.gammaln(s + 1) - sp.gammaln(x + 1)) + mp.exp(
					np.log(1. - w) + np.log(2) + (s + x) * np.log(L) + x * np.log(r) + (0.5 * (-s - x + a)) * np.log(
						(L * (1 + r)) / b) + a * np.log(b) + mp.log(
						mp.besselk(s + x - a, 2 * np.sqrt(L * (1 + r) * b))) - sp.gammaln(s + 1) - sp.gammaln(
						x + 1) - sp.gammaln(a))) / pofs(s, L, thr)
			else:
				return (np.exp(np.log(w) + (s + x) * np.log(L) + x * np.log(r) + np.log(t) + (-1 - s - x) * np.log(
					L + L * r + t) + sp.gammaln(1 + s + x) - sp.gammaln(s + 1) - sp.gammaln(x + 1)) + np.exp(
					np.log(1. - w) + np.log(2) + (s + x) * np.log(L) + x * np.log(r) + (0.5 * (-s - x + a)) * np.log(
						(L * (1 + r)) / b) + a * np.log(b) + np.log(
						sp.kv(s + x - a, 2 * np.sqrt(L * (1 + r) * b))) - sp.gammaln(s + 1) - sp.gammaln(
						x + 1) - sp.gammaln(a))) / pofs(s, L, thr)
	elif modC == 5:
		# *************** lambda ~ w * Gamma + (1-w) * Gamma (Gamma mixture model):
		def pofs(s, L):
			return np.exp(np.log(w) + s * np.log(L * b) + (-s - a) * np.log(1 + L * b) + sp.gammaln(s + a) - sp.gammaln(
				s + 1) - sp.gammaln(a)) + np.exp(
				np.log(1. - w) + s * np.log(L * d) + (-s - g) * np.log(1 + L * d) + sp.gammaln(s + g) - sp.gammaln(
					s + 1) - sp.gammaln(g))

		def pofx_given_s(x, s, L, r, thr):
			return (np.exp(np.log(w) + x * np.log(r) + (s + x) * np.log(L * b) + (-s - x - a) * np.log(
				1 + L * (1 + r) * b) + sp.gammaln(s + x + a) - sp.gammaln(s + 1) - sp.gammaln(x + 1) - sp.gammaln(
				a)) + np.exp(np.log(1 - w) + x * np.log(r) + (s + x) * np.log(L * d) + (-s - x - g) * np.log(
				1 + L * (1 + r) * d) + sp.gammaln(s + x + g) - sp.gammaln(s + 1) - sp.gammaln(x + 1) - sp.gammaln(
				g))) / pofs(s, L)
	elif modC == 6:
		# *************** lambda ~ w * Gamma + (1-w) * InvGamma (mixture model):
		def pofs(s, L, thr):
			if thr:
				return np.exp(
					np.log(w) + s * np.log(L * b) + (-s - a) * np.log(1 + L * b) + sp.gammaln(s + a) - sp.gammaln(
						s + 1) - sp.gammaln(a)) + mp.exp(
					np.log(1. - w) + np.log(2.) + ((s + g) / 2.) * np.log(L * d) + mp.log(
						mp.besselk(-s + g, 2 * mp.sqrt(L * d))) - sp.gammaln(s + 1) - sp.gammaln(g))
			else:
				return np.exp(
					np.log(w) + s * np.log(L * b) + (-s - a) * np.log(1 + L * b) + sp.gammaln(s + a) - sp.gammaln(
						s + 1) - sp.gammaln(a)) + np.exp(
					np.log(1. - w) + np.log(2.) + ((s + g) / 2.) * np.log(L * d) + np.log(
						sp.kv(-s + g, 2 * np.sqrt(L * d))) - sp.gammaln(s + 1) - sp.gammaln(g))

		def pofx_given_s(x, s, L, r, thr):
			if thr:
				return (np.exp(np.log(w) + x * np.log(r) + (s + x) * np.log(L * b) + (-s - x - a) * np.log(
					1 + L * (1 + r) * b) + sp.gammaln(s + x + a) - sp.gammaln(s + 1) - sp.gammaln(x + 1) - sp.gammaln(
					a)) + mp.exp(
					np.log(1 - w) + np.log(2) + (s + x) * np.log(L) + x * np.log(r) + (0.5 * (-s - x + g)) * np.log(
						(L * (1 + r)) / d) + g * np.log(d) + mp.log(
						mp.besselk(s + x - g, 2 * np.sqrt(L * (1 + r) * d))) - sp.gammaln(s + 1) - sp.gammaln(
						x + 1) - sp.gammaln(g))) / pofs(s, L, thr)
			else:
				return (np.exp(np.log(w) + x * np.log(r) + (s + x) * np.log(L * b) + (-s - x - a) * np.log(
					1 + L * (1 + r) * b) + sp.gammaln(s + x + a) - sp.gammaln(s + 1) - sp.gammaln(x + 1) - sp.gammaln(
					a)) + np.exp(
					np.log(1 - w) + np.log(2) + (s + x) * np.log(L) + x * np.log(r) + (0.5 * (-s - x + g)) * np.log(
						(L * (1 + r)) / d) + g * np.log(d) + np.log(
						sp.kv(s + x - g, 2 * np.sqrt(L * (1 + r) * d))) - sp.gammaln(s + 1) - sp.gammaln(
						x + 1) - sp.gammaln(g))) / pofs(s, L, thr)

	pvals = []
	L = 1.
	gcnt = 0
	for gene in genes:
		gcnt += 1
		if gcnt % 1000 == 0:
			sys.stderr.write("%i%% done.\r" % int(float(gcnt)/len(genes)*100.))

		sobs = int(gene["obs"][2])
		mobs = int(gene["obs"][0])
		kobs = int(gene["obs"][1])
		sexp = gene["exp"][2]
		mexp = gene["exp"][0]
		kexp = gene["exp"][1]
		ratm = mexp / sexp
		ratk = kexp / sexp

		large_flag = 0
		last_p = 0.
		if sobs == 0:
			meant2 = 2  # ~= E[m] * 2
		else:
			meant2 = int(ratm * sobs + 3. * math.sqrt(ratm * sobs))	# ~= E[m] + 3*sigma_m
		for mtest in range(meant2):
			if math.isnan(pofx_given_s(mtest, sobs, L, ratm, 0)):
				large_flag = 1  # Going to large-x mode.

				class P_of_x_given_s(st.rv_continuous):
					def _pdf(self, x, s, eL, rat):
						if s == 1e-10:
							s = 0
						return pofx_given_s(x, s, eL, rat, 1)

				inst_pofx = P_of_x_given_s(a=0)
				break
			cur_p = pofx_given_s(mtest, sobs, L, ratm, 0)
			diff = cur_p - last_p
			last_p = cur_p
			if pofx_given_s(mtest, sobs, L, ratm, 0) > 1.:
				if pofx_given_s(mtest - 1, sobs, L, ratm, 0) > 1. / runC or diff > 0:
					large_flag = 1  # Going to large-x mode.

					class P_of_x_given_s(st.rv_continuous):
						def _pdf(self, x, s, eL, rat):
							if s == 1e-10:
								s = 0
							return pofx_given_s(x, s, eL, rat, 1)

					inst_pofx = P_of_x_given_s(a=0)
				else:
					class P_of_x_given_s(st.rv_discrete):
						def _pmf(self, x, s, eL, rat):
							if s == 1e-10:
								s = 0
							return pofx_given_s(x, s, eL, rat, 0)

					inst_pofx = P_of_x_given_s()
				break
		if large_flag == 0:
			class P_of_x_given_s(st.rv_discrete):
				def _pmf(self, x, s, eL, rat):
					if s == 1e-10:
						s = 0
					return pofx_given_s(x, s, eL, rat, 0)

			inst_pofx = P_of_x_given_s()

		# The upper limit for mtest should coincide with the upper limit above (meant2).
		sum_p = 0.
		testm_array = []
		for mtest in range(meant2):
			m_pneg = sum_p + pofx_given_s(mtest, sobs, L, ratm, large_flag)
			m_ppos = 1. - sum_p
			m_ppos = m_ppos.real
			m_pneg = m_pneg.real
			if m_ppos < 0 or m_ppos > 1:
				sys.stderr.write("m_ppos outside limits: %i \t %e. \n" % (large_flag, m_ppos)) 
			testm_array.append([m_pneg, m_ppos])
			sum_p += pofx_given_s(mtest, sobs, L, ratm, large_flag)

		sum_p = 0.
		testk_array = []
		for ktest in range((int(ratk * sobs) + 1) * 2):
			k_pneg = sum_p + pofx_given_s(ktest, sobs, L, ratk, large_flag)
			k_ppos = 1. - sum_p
			k_ppos = k_ppos.real
			k_pneg = k_pneg.real
			if k_ppos < 0 or k_ppos > 1:
				sys.stderr.write("k_ppos outside limits: %i \t %e. \n" % (large_flag, k_ppos)) 
			testk_array.append([k_pneg, k_ppos])
			sum_p += pofx_given_s(ktest, sobs, L, ratk, large_flag)

		for rep in range(runC):
			if simC:
				#	Simulate expectation under null.
				if sobs == 0:
					sobs = 1e-10
				try:
					msim = inst_pofx.rvs(sobs, L, ratm)
					ksim = inst_pofx.rvs(sobs, L, ratk)
				except:
					continue  # If this happens *very* frequently (out of runC*18,666 runs in total), try decreasing meant2 above. Otherwise results may overestimate the degree of negative selection.
				if sobs == 1e-10:
					sobs = 0
				mobs = int(round(msim))
				kobs = int(round(ksim))

			try:
				[m_pneg, m_ppos] = testm_array[mobs]
				[k_pneg, k_ppos] = testk_array[kobs]
			except:
				cum_p = 0.
				for x in range(mobs):
					cum_p += pofx_given_s(x, sobs, L, ratm, large_flag)
				m_pneg = cum_p + pofx_given_s(mobs, sobs, L, ratm, large_flag)
				m_ppos = 1. - cum_p
				m_ppos = float(m_ppos.real)
				m_pneg = float(m_pneg.real)
				if m_ppos < 0 or m_ppos > 1 or math.isinf(m_ppos) or math.isnan(m_ppos):
					cum_p = 0.
					for x in range(mobs):
						cum_p += pofx_given_s(x, sobs, L, ratm, 1)
					m_pneg = cum_p + pofx_given_s(mobs, sobs, L, ratm, 1)
					m_ppos = 1. - cum_p
					m_ppos = float(m_ppos.real)
					m_pneg = float(m_pneg.real)
					if m_ppos < 0 or m_ppos > 1 or math.isinf(m_ppos) or math.isnan(m_ppos):
						sys.stderr.write("Setting p_m^pos --> 0 on gene %s (was %e).\n" % (gene["gene"], m_ppos))
						m_ppos = 0
				if m_pneg > 1:
					if m_pneg < 1.01:
						m_pneg = 1.
					elif math.isinf(m_pneg) or math.isnan(m_pneg):
						if simC == 0:
							sys.stderr.write(
								"p_m^neg on gene %s: %f --> 1.\n" % (gene["gene"], m_pneg))
						m_pneg = 1.

				cum_p = 0.
				for x in range(kobs):
					cum_p += pofx_given_s(x, sobs, L, ratk, large_flag)
				k_pneg = cum_p + pofx_given_s(kobs, sobs, L, ratk, large_flag)
				k_ppos = 1. - cum_p
				k_ppos = float(k_ppos.real)
				k_pneg = float(k_pneg.real)
				if k_ppos < 0 or k_ppos > 1 or math.isinf(k_ppos) or math.isnan(k_ppos):
					cum_p = 0.
					for x in range(kobs):
						cum_p += pofx_given_s(x, sobs, L, ratk, 1)
					k_pneg = cum_p + pofx_given_s(kobs, sobs, L, ratk, 1)
					k_ppos = 1. - cum_p
					k_ppos = float(k_ppos.real)
					k_pneg = float(k_pneg.real)
					if k_ppos < 0 or k_ppos > 1 or math.isinf(k_ppos) or math.isnan(k_ppos):
						sys.stderr.write("Setting p_k^pos --> 0 on gene %s (was %e).\n" % (gene["gene"], k_ppos))
						k_ppos = 0
				if k_pneg > 1:
					if k_pneg < 1.01:
						k_pneg = 1.
					elif math.isinf(k_pneg) or math.isnan(k_pneg):
						if simC == 0:
							sys.stderr.write(
								"p_k^neg on gene %s: %f --> 1.\n" % (gene["gene"], k_pneg))
						k_pneg = 1.

			pvals.append([gene["gene"], m_pneg, k_pneg, m_ppos, k_ppos, pofx_given_s(0, sobs, L, ratm, 0),
						  pofx_given_s(0, sobs, L, ratk, 0), mobs, kobs, sobs])

	sys.stderr.write("100% done.\n")

	return pvals


def construct_histogram(var_array, bin_var):
	noinf = [el for el in var_array if el < 1e50]
	var_max = max(noinf) + bin_var
	hist = [0. for i in range(int(var_max/bin_var))]
	for var in noinf:
		hist[int(var / bin_var)] += 1
	return [[(i + 0.5) * bin_var, hist[i] / len(noinf)] for i in range(len(hist))]


def compute_phi_sim(pvals_array, ind1, ind2):
	all_phi = []
	for gene in pvals_array:
		if gene[ind1] == 0. or gene[ind2] == 0.:
			all_phi.append(1e5)
		else:
			all_phi.append(-math.log(gene[ind1]) - math.log(gene[ind2]))
	return all_phi


def compute_phi_obs(pvals_array, ind1, ind2):
	all_phi = []
	for gene in pvals_array:
		if gene[ind1] == 0. or gene[ind2] == 0.:
			cur_phi = 1e5
		else:
			cur_phi = -math.log(gene[ind1]) - math.log(gene[ind2])
		all_phi.append(
			{"gene": gene[0], "phi": cur_phi, "p0m": gene[5], "p0k": gene[6], "mks": [gene[7], gene[8], gene[9]]})
	return all_phi


def fdr_bh(pvals_obs):
	raise NotImplementedError


def FDR_discrete(phi_sim_array, gene_phi_real, bin_phi, bin_p):
	phi_sim_hist = construct_histogram(phi_sim_array, bin_phi)

	gene_phi_real.sort(key=lambda arg: arg["phi"], reverse=True)

	larger_0 = [el for el in gene_phi_real if el["phi"] > 0.]
	equal_0 = sorted([el for el in gene_phi_real if abs(el["phi"]) < 1e-30], key=lambda arg: arg["p0m"] * arg["p0k"],
					 reverse=True)
	sorted_phi = larger_0 + equal_0
	if len(sorted_phi) < len(gene_phi_real):
		excl = [el for el in gene_phi_real if [gen["gene"]
											   for gen in sorted_phi].count(el["gene"]) == 0]
		sys.stderr.write("Warning: Not including gene(s) in output:\n")
		for el in excl:
			sys.stderr.write("%s\n" % el)

	phi_sim_max = max(phi_sim_array)

	phi_pvals_obs = []
	for gene in sorted_phi:
		if gene["phi"] > phi_sim_max + 0.5 * bin_phi:
			phi_pvals_obs.append({"gene": gene["gene"], "p_phi": 0, "phi": gene["phi"], "mks": gene["mks"]})
		else:
			if abs(gene["phi"]) < 1e-30:
				phi_pvals_obs.append(
					{"gene": gene["gene"], "p_phi": 1., "phi": gene["p0m"] * gene["p0k"], "mks": gene["mks"]})
			else:
				cumprob = 0.
				i = 0
				while i < len(phi_sim_hist) and phi_sim_hist[i][0] + 0.5 * bin_phi <= max(0., gene["phi"]):
					cumprob += phi_sim_hist[i][1]
					i += 1
				phi_pvals_obs.append(
					{"gene": gene["gene"], "p_phi": 1. - cumprob, "phi": gene["phi"], "mks": gene["mks"]})

	phi_pvals_sim = []
	for phis in phi_sim_array:
		cumprob = 0.
		i = 0
		while i < len(phi_sim_hist) and phi_sim_hist[i][0] + 0.5 * bin_phi <= max(0., phis):
			cumprob += phi_sim_hist[i][1]
			i += 1
		phi_pvals_sim.append(1. - cumprob)

	p_sim_hist = construct_histogram(phi_pvals_sim, bin_p)

	phi_qvals_obs = []
	cumprob = 0.
	i = 0
	gcnt = 0.
	for gene in phi_pvals_obs:
		gcnt += 1
		if gene["p_phi"] == 0:
			phi_qvals_obs.append(
				{"gene": gene["gene"], "p_phi": gene["p_phi"], "q_phi": 0, "phi": gene["phi"], "mks": gene["mks"]})
		else:
			while i < len(p_sim_hist) and p_sim_hist[i][0] - 0.5 * bin_p <= max(0., gene["p_phi"]):
				cumprob += p_sim_hist[i][1]
				i += 1
			phi_qvals_obs.append(
				{"gene": gene["gene"], "p_phi": gene["p_phi"], "q_phi": cumprob / gcnt * len(phi_pvals_obs),
				 "phi": gene["phi"], "mks": gene["mks"]})

	if len(phi_qvals_obs) != len(gene_phi_real):
		sys.stderr.write("Number of genes in output different from original: %i vs. %i.\n" % (
			len(phi_qvals_obs), len(gene_phi_real)))

	phi_qvals_adj = []
	for g in range(len(phi_qvals_obs)):
		cur_min = min([el["q_phi"] for el in phi_qvals_obs[g:]])
		gene = phi_qvals_obs[g]
		phi_qvals_adj.append(
			{"gene": gene["gene"], "p_phi": gene["p_phi"], "q_adj": min(1., cur_min), "phi": gene["phi"],
			 "mks": gene["mks"]})

	return phi_qvals_adj

# ************************************* GLOBAL DEFINITIONS **************************************


mod_choice = ["", "Gamma(a,b): [a,b] =", "InverseGamma(a,b): [a,b] =", "w Exp(t) + (1-w) Gamma(a,b): [a,b,t,w] =",
			  "w Exp(t) + (1-w) InverseGamma(a,b): [a,b,t,w] =", "w Gamma(a,b) + (1-w) Gamma(g,d): [a,b,g,d,w] =",
			  "w Gamma(a,b) + (1-w) InverseGamma(g,d): [a,b,g,d,w] ="]  # model choice
modC_map = [2, 2, 4, 4, 5, 5]			# map model --> number of params
run_no = 100  # No. of simulation replicates used for computing FDR
outname = str(sys.argv[1])  # name for the output file containing the q values

# ***********************************************************************************************
# ***********************************************************************************************
#
# Collect parameter estimates from all fitted models in working directory.

p_files = glob.glob("Output/param_estimates_%s_*.txt" % outname)

repeat = True
while repeat: 
	try:
		all_models = []
		for p_file in p_files:
			fin = open(p_file)
			lines = fin.readlines()
			fin.close()
			field = lines[0].strip().split(", ")
			all_models.append([float(el) for el in field])

		cur_min = 1e20
		cur_ind = 10
		for m in range(len(all_models)):
			if 2.*modC_map[int(all_models[m][-1])-1] + 2.*all_models[m][-2] < cur_min:
				cur_min = 2.*modC_map[int(all_models[m][-1])-1] + 2.*all_models[m][-2]
				cur_ind = m
		if cur_min < 1e20:
			sys.stderr.write("Best model fit: model %i.\n" %
							int(all_models[cur_ind][-1]))
			fout = open("Output/used_params_and_model_%s.txt" % outname, "w")
			fout.write(''.join([''.join(["%e, " for i in range(modC_map[int(all_models[cur_ind][-1])-1])]),
					"%i\n"]) % tuple(all_models[cur_ind][:-2] + [int(all_models[cur_ind][-1])]))
			fout.close()
		else:
			sys.stderr.write("Could not find a converging solution.\n")


		# ***********************************************************************************************
		#
		# Compute q-values for all genes, using best fitting model.

		# Import parameters and index of chosen model.
		fin = open("Output/used_params_and_model_%s.txt" % outname)
		lines = fin.readlines()
		fin.close()
		field = lines[0].strip().split(", ")

		mod_C = int(field[-1])
		params = [float(el) for el in field[:-1]]

		if len(params) != modC_map[mod_C-1]:
			sys.stderr.write("Number of inferred parameters does not match the chosen model: %i vs. %i.\n" % (
				len(params), modC_map[mod_C-1]))
			sys.exit()

		# Import output from data_preparation, containing l_x and (m,k,s)_obs.
		fin = open("Output/output_data_preparation_%s.txt" % outname)
		lines = fin.readlines()
		fin.close()
		mks_type = []
		for line in lines:
			field = line.strip().split("\t")
			mks_type.append({"gene": field[0], "exp": [float(field[1]), float(field[2]), float(field[3])],
							"obs": [float(field[4]), float(field[5]), float(field[6])], "len": int(field[7])})
		mks_type = sorted(mks_type, key=lambda arg: arg["gene"])

		sys.stderr.write("Computing real p-values.\n")
		pvals_obs = compute_p_values(params, mks_type, [mod_C, 0, 1])
		sys.stderr.write("Computing simulated p-values.\n")
		sys.stderr.write("Simulation runs\t= %i.\n" % (run_no))
		pvals_sim = compute_p_values(params, mks_type, [mod_C, 1, run_no])
		# Format: [gene, m_pneg, k_pneg, m_ppos, k_ppos, pofx_given_s(0, sobs, L, ratm), pofx_given_s(0, sobs, L, ratk)]

		#	Negative selection
		# phi_sim = compute_phi_sim(pvals_sim, 1, 2)
		# gene_phi_obs = compute_phi_obs(pvals_obs, 1, 2)
		# q_neg_adj = FDR_discrete(phi_sim, gene_phi_obs, 0.02, 0.000001)

		#	Positive selection
		phi_sim = compute_phi_sim(pvals_sim, 3, 4)
		gene_phi_obs = compute_phi_obs(pvals_obs, 3, 4)
		q_pos_adj = FDR_discrete(phi_sim, gene_phi_obs, 0.02, 0.000001)

		# Sort and output q-values.
		# q_neg_sorted = sorted(q_neg_adj, key=lambda arg: arg["gene"])
		q_pos_sorted = sorted(q_pos_adj, key=lambda arg: arg["gene"])
		pvals = [val["p_phi"] for val in q_pos_sorted]
		qvals = multipletests(pvals, method='fdr_bh', )[1]

		repeat = False
	except Exception as e:
		sys.stderr.write(str(e))
		
		repeat = False

fout = open("Output/q_values_%s.txt" % outname, "w")
fout.write("%s\t%s\n" % (mod_choice[mod_C], params))
fout.write(
	"gene\tp_pos\tq_pos\tmis_obs\tnon_obs\tsyn_obs\n")
for g in range(len(q_pos_sorted)):
	fout.write("%s\t%e\t%e\t%i\t%i\t%i\n"
			   % (q_pos_sorted[g]["gene"], pvals[g], qvals[g], q_pos_sorted[g]["mks"][0],
				  q_pos_sorted[g]["mks"][1], q_pos_sorted[g]["mks"][2]))
fout.close()

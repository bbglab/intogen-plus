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
# * Version:		1.1										 *
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
import itertools as it
import math
import gzip
import numpy as np
import random
import mpmath as mp
import scipy.stats as st
import itertools as it


# ************************************ FUNCTION DEFINITIONS *************************************

def muttype_index(cod_ID):
	# 0: missense, 1: nonsense, 2: synonymous
	if [56, 57, 58, 59].count(cod_ID):
		return [[0, 0, 0, 0], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif [20, 21, 22, 23].count(cod_ID):
		return [[0, 0, 0, 0], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif [16, 17, 18, 19].count(cod_ID):
		return [[0, 0, 0, 0], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif [24, 25, 26, 27].count(cod_ID):
		return [[0, 0, 0, 0], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif cod_ID == 40:
		return [[0, 0, 0, 1], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif [41, 42, 43].count(cod_ID):
		return [[0, 0, 0, 0], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif cod_ID == 5:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 2]]
	elif cod_ID == 7:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0]]
	elif cod_ID == 1:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 2]]
	elif cod_ID == 3:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0]]
	elif cod_ID == 9:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 2]]
	elif cod_ID == 11:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0]]
	elif cod_ID == 61:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 2]]
	elif cod_ID == 63:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0]]
	elif cod_ID == 13:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [1, 0, 1, 2]]
	elif cod_ID == 15:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [1, 2, 1, 0]]
	elif cod_ID == 4:
		return [[0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 2, 0]]
	elif cod_ID == 6:
		return [[0, 0, 0, 1], [0, 0, 0, 0], [2, 0, 0, 0]]
	elif cod_ID == 0:
		return [[0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 2, 0]]
	elif cod_ID == 2:
		return [[0, 0, 0, 1], [0, 0, 0, 0], [2, 0, 0, 0]]
	elif cod_ID == 8:
		return [[0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 2, 0]]
	elif cod_ID == 10:
		return [[0, 0, 0, 1], [0, 0, 0, 0], [2, 0, 0, 0]]
	elif cod_ID == 45:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 2]]
	elif cod_ID == 47:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [1, 2, 0, 0]]
	elif [52, 54].count(cod_ID):
		return [[0, 0, 0, 2], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif [53, 55].count(cod_ID):
		return [[0, 0, 0, 0], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif cod_ID == 60:
		return [[0, 2, 0, 0], [1, 0, 1, 0], [0, 0, 2, 0]]
	elif cod_ID == 62:
		return [[0, 2, 0, 0], [1, 0, 0, 0], [2, 0, 0, 0]]
	elif cod_ID == 36:
		return [[2, 0, 0, 1], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif cod_ID == 38:
		return [[2, 0, 0, 0], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif [37, 39].count(cod_ID):
		return [[0, 0, 0, 0], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif cod_ID == 32:
		return [[0, 2, 0, 1], [0, 0, 0, 0], [0, 0, 2, 0]]
	elif cod_ID == 34:
		return [[0, 2, 0, 0], [0, 0, 0, 0], [2, 0, 0, 0]]
	elif cod_ID == 28:
		return [[0, 0, 0, 0], [1, 0, 1, 0], [2, 2, 2, 2]]
	elif cod_ID == 30:
		return [[0, 0, 0, 0], [1, 0, 0, 0], [2, 2, 2, 2]]
	elif [29, 31].count(cod_ID):
		return [[0, 0, 0, 0], [0, 0, 0, 0], [2, 2, 2, 2]]
	elif cod_ID == 33:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 2]]
	elif cod_ID == 35:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0]]
	elif cod_ID == 12:
		return [[1, 1, 1, 1], [0, 1, 2, 1], [0, 1, 2, 1]]
	elif cod_ID == 14:
		return [[1, 1, 1, 1], [1, 1, 1, 1], [2, 1, 0, 1]]
	elif cod_ID == 44:
		return [[1, 1, 1, 1], [2, 1, 0, 1], [1, 1, 1, 1]]
	elif [48, 49, 51].count(cod_ID):
		return [[0, 0, 0, 0], [0, 0, 0, 0], [2, 2, 0, 2]]
	elif cod_ID == 50:
		return [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
	elif cod_ID == 46:
		return [[0, 0, 0, 0], [1, 0, 0, 0], [1, 0, 0, 0]]
	else:
		sys.stderr.write("Strange codon ID: %i.\n" % cod_ID)


def import_special_genes(filename):
	fin = open(filename)
	lines = fin.readlines()
	fin.close()
	c_genes = []
	for line in lines:
		c_genes.append(line.strip().split()[0])
	return c_genes


def import_context_freqs(filename):
	fin = open(filename)
	lines = fin.readlines()
	fin.close()
	occs = []
	for line in lines:
		field = line.strip().split()
		occs.append(float(field[1]))
	return occs


def import_quintuplets(filename):
	fin = open(filename)
	lines = fin.readlines()
	fin.close()
	return [line.strip().split()[0] for line in lines]


def import_maf_data(filename, context_mode):
	fin = open(filename)
	lines = fin.readlines()
	fin.close()
	mut_array = []
	for line in lines[1:]:
		field = line.strip().split("\t")
		if len(field) != 4:
			sys.stderr.write("Number of columns in maf file not as expected (=4): %i.\n" % len(field))
			sys.exit()
		# context is 0-based; triplets mapped to legacy
		if context_mode == 0:
			mut_array.append({"gene": field[0], "muttype": field[1], "mutbase": field[2].upper(),
							  "context": triplets.index(triplets_user[int(field[3])])})
		else:
			mut_array.append(
				{"gene": field[0], "muttype": field[1], "mutbase": field[2].upper(), "context": int(field[3])})
	return mut_array


def import_known_genes_UCSC(filename):
	fin = open(filename)
	lines = fin.readlines()
	fin.close()
	genes = []
	for line in lines:
		field = line.strip().split('\t')
		if field[1][3:] == 'M' or len(field[1].split('_')) > 1:
			continue
		if field[1][3:] == 'X':
			chrnr = 23
		elif field[1][3:] == 'Y':
			chrnr = 24
		else:
			chrnr = int(field[1][3:])
		if len(field) != 14:
			sys.stderr.write("Unexpected line format.\n")
			sys.exit()
		if field[11] != 'n/a' and field[12] != 'n/a':
			genes.append({"transcript": field[0], "gene": field[10], "chr": chrnr, "strand": field[2],
						  "genebegin": int(field[3]), "geneend": int(field[4]), "cdbegin": int(field[5]),
						  "cdend": int(field[6]), "exoncnt": int(field[7]),
						  "exonbegins": [int(el) for el in field[8][:-1].split(',')],
						  "exonends": [int(el) for el in field[9][:-1].split(',')], "pepseq": field[12],
						  "ensemblID": field[13], "all": field})
	genes = sorted(genes, key=lambda arg: arg["gene"])
	return genes


def import_codons_by_gene(filename):
	sys.stderr.write("Importing codons...\n")
	with gzip.open(filename, 'rt') as fin:
		lines = fin.readlines()
	fin.close()
	c_genes = []
	cur_gene = "bla"
	cur_cods = []
	for line in lines:
		field = line.strip().split("\t")
		if len(field) == 2:
			c_genes.append({"gene": cur_gene, "context_triplets": cur_cods})
			cur_cods = []
			cur_gene = field[1]
		else:
			cur_cods.append([int(el) for el in field])
	c_genes.append({"gene": cur_gene, "context_triplets": cur_cods})
	return c_genes[1:]


def make_neutral_mut_matrix_pentanucs(neut_cont_array, penta_occs_array):
	neutral_mut_matrix = [[0. for i in range(4)] for j in range(1024)]
	totcnt = 0
	for el in neut_cont_array:
		neutral_mut_matrix[el[0]][el[1]] += 1.
		totcnt += 1
	#	Including pseudocounts:
	for m in range(16):
		for k in range(4):
			alist = [0, 1, 2, 3]
			alist.remove(k)
			uselist = alist[:]
			for i in range(16):
				ind1 = m * 64 + k * 16 + i
				for j in uselist:
					neutral_mut_matrix[ind1][j] += 1.
					totcnt += 1
	# Symmetrize the probabilities:
	probs_array = [[0. for i in range(4)] for j in range(1024)]
	for i in range(4):
		for j in range(4):
			for k in range(4):
				for l in range(4):
					for m in range(4):
						n = 256 * i + 64 * j + 16 * k + 4 * l + m
						nprime = 1023 - (256 * m + 64 * l + 16 * k + 4 * j + i)
						for h in range(4):
							totalcnt = neutral_mut_matrix[n][h] + neutral_mut_matrix[nprime][3 - h]
							probs_array[n][h] = totalcnt / (2. * totcnt) / (
										(penta_occs_array[n] + penta_occs_array[nprime]) / 2.)
							probs_array[nprime][3 - h] = probs_array[n][h]
	return probs_array


def make_neutral_mut_matrix_trinucs(neut_cont_array, trinuc_occs_array):
	neutral_mut_matrix = [[0. for i in range(4)] for j in range(64)]
	totcnt = 0
	for el in neut_cont_array:
		neutral_mut_matrix[el[0]][el[1]] += 1.
		totcnt += 1
	# Symmetrize the probabilities:
	probs_array = [[0. for i in range(4)] for j in range(64)]
	for i in range(4):
		for j in range(4):
			for k in range(4):
				n = 16 * j + 4 * i + k
				nprime = 63 - (16 * j + 4 * k + i)
				for m in range(4):
					totalcnt = neutral_mut_matrix[n][m] + neutral_mut_matrix[nprime][3 - m]
					probs_array[n][m] = totalcnt / (2. * totcnt) / (
								(trinuc_occs_array[n] + trinuc_occs_array[nprime]) / 2.)
					probs_array[nprime][3 - m] = probs_array[n][m]
	return probs_array


def export_expected_observed_mks_per_gene(codon_array_by_gene, mut_array, neutral_muts_array, nuc_occs_array,
										  context_mode):
	mut_coding = [m for m in mut_array if ["missense", "nonsense", "coding-synon"].count(m["muttype"])]
	sys.stderr.write("%i (m,k,s) mutations.\n" % (len(mut_coding)))

	gene_keys = []
	muts_by_gene = []
	for k, g in it.groupby(sorted(mut_coding, key=lambda arg: arg["gene"]), key=lambda arg: arg["gene"]):
		gene_keys.append(k)
		muts_by_gene.append(list(g))

	#	Derive the neutral transition matrix from the sum over all patients.
	if context_mode == 0:  # TRINUCLEOTIDES
		neutral_matrix = make_neutral_mut_matrix_trinucs(neutral_muts_array, nuc_occs_array)
	else:  # PENTANUCLEOTIDES
		neutral_matrix = make_neutral_mut_matrix_pentanucs(neutral_muts_array, nuc_occs_array)
	for el in neutral_matrix:
		for i in range(len(el)):
			sys.stderr.write("%f\t" % el[i])
		sys.stderr.write("\n")
	if len(neutral_matrix) == 0:
		sys.stderr.write("Cannot construct neutral mutation matrix.\n")
		sys.exit()

	exp_obs_per_gene = []
	for g in range(len(codon_array_by_gene)):

		cur_gene = codon_array_by_gene[g]["gene"]
		codons_gene = codon_array_by_gene[g]["context_triplets"]
		gene_len = len(codons_gene) * 3

		xobs = [0 for i in range(3)]
		try:
			gene_muts = muts_by_gene[gene_keys.index(cur_gene)]
			#	Sum observed mutations in categories over patients.
			for mut in gene_muts:
				xobs[["missense", "nonsense", "coding-synon"].index(mut["muttype"])] += 1
		except:
			pass

		#	Compute neutral expectation (up to a constant factor).
		expect_x = [0. for t in range(3)]
		cod_cnt = 0
		for codon in codons_gene:
			for i in range(3):
				if context_mode:  # PENTANUCLEOTIDES
					if i == 0:
						if cod_cnt == 0:
							continue
						else:
							cur_pent = triplets[codons_gene[cod_cnt - 1][2]] + triplets[codon[2]][:2]
					elif i == 1:
						cur_pent = triplets[codon[0]] + triplets[codon[2]][1:]
					else:
						if cod_cnt == len(codons_gene) - 1:
							continue
						else:
							cur_pent = triplets[codon[1]] + triplets[codons_gene[cod_cnt + 1][0]][1:]
				for k in range(4):
					codon_ID = codon[1]
					if context_mode == 0:  # TRINUCLEOTIDES
						expect_x[muttype_index(codon_ID)[i][k]] += neutral_matrix[codon[i]][k]
					else:  # PENTANUCLEOTIDES
						expect_x[muttype_index(codon_ID)[i][k]] += neutral_matrix[quintuplets.index(cur_pent)][k]
			cod_cnt += 1
		exp_obs_per_gene.append([cur_gene, expect_x[0], expect_x[1], expect_x[2], xobs[0], xobs[1], xobs[2], gene_len])
	return exp_obs_per_gene


def minimize_neg_ln_L(p_start, function, mks_array, aux, bound_array, n_param):
	if n_param == 2:
		p0, p1 = p_start
		res = minimize(function, (p0, p1), args=(mks_array, aux), method='L-BFGS-B', bounds=bound_array,
					   options={'disp': None, 'gtol': 1e-12, 'eps': 1e-5, 'maxiter': 15000, 'ftol': 1e-12})
		return [res.x[0], res.x[1], res.fun]
	elif n_param == 4:
		p0, p1, p2, p3 = p_start
		res = minimize(function, (p0, p1, p2, p3), args=(mks_array, aux), method='L-BFGS-B', bounds=bound_array,
					   options={'disp': None, 'gtol': 1e-12, 'eps': 1e-5, 'maxiter': 15000, 'ftol': 1e-12})
		return [res.x[0], res.x[1], res.x[2], res.x[3], res.fun]
	elif n_param == 5:
		p0, p1, p2, p3, p4 = p_start
		res = minimize(function, (p0, p1, p2, p3, p4), args=(mks_array, aux), method='L-BFGS-B', bounds=bound_array,
					   options={'disp': None, 'gtol': 1e-12, 'eps': 1e-5, 'maxiter': 15000, 'ftol': 1e-12})
		return [res.x[0], res.x[1], res.x[2], res.x[3], res.x[4], res.fun]


def neg_ln_L(p, genes, aux):
	modC = aux
	if [1, 2].count(modC):
		a, b = p
		if a < 0:
			a = 1e-6
		if b < 0:
			b = 1e-6
	elif [3, 4].count(modC):
		a, b, t, w = p
		if a < 0:
			a = 1e-6
		if b < 0:
			b = 1e-6
		if w < 0:
			w = 1e-6
		if t < 0:
			t = 1e-6
	elif [5, 6].count(modC):
		a, b, g, d, w = p
		if a < 0:
			a = 1e-6
		if b < 0:
			b = 1e-6
		if g < 0:
			g = 1e-6
		if d < 0:
			d = 1e-6
		if w < 0:
			w = 1e-6

	genes_by_sobs = [[ka, len(list(gr))] for ka, gr in it.groupby(sorted(
		genes, key=lambda arg: int(arg["obs"][2])), key=lambda arg: int(arg["obs"][2]))]

	summe = 0.
	if modC == 1:
		for sval in genes_by_sobs:
			s = sval[0]
			# *************** lambda ~ Gamma:
			summe += sval[1] * (s * np.log(b) + (-s - a) * np.log(1 + b) + sp.gammaln(s + a) - sp.gammaln(s + 1) - sp.gammaln(a))
	elif modC == 2:
		for sval in genes_by_sobs:
			s = sval[0]
			# *************** lambda ~ IG:
			if s > 25:
				summe += sval[1] * (math.log(2.) + ((s + a) / 2.) * math.log(b) + float(mp.log(
					mp.besselk(-s + a, 2 * math.sqrt(b)).real)) - sp.gammaln(s + 1) - sp.gammaln(a))
			else:
				try:
					summe += sval[1] * (math.log(2.) + ((s + a) / 2.) * math.log(b) + math.log(
						sp.kv(-s + a, 2 * math.sqrt(b))) - sp.gammaln(s + 1) - sp.gammaln(a))
				except:
					summe += sval[1] * (math.log(2.) + ((s + a)/2.)*math.log(b) + float(mp.log(
						mp.besselk(-s + a, 2 * math.sqrt(b)).real)) - sp.gammaln(s + 1) - sp.gammaln(a))
	elif modC == 3:
		for sval in genes_by_sobs:
			s = sval[0]
			# *************** lambda ~ w * Exp + (1-w) * Gamma:
			summe += sval[1] * (np.log(math.exp(math.log(w * t) + (-1 - s) * math.log(1 + t)) + math.exp(
				math.log(1. - w) + s * math.log(b) + (-s - a) * math.log(1 + b) + sp.gammaln(s + a) - sp.gammaln(
					s + 1) - sp.gammaln(a))))
	elif modC == 4:
		for sval in genes_by_sobs:
			s = sval[0]
			# *************** lambda ~ w * Exp + (1-w) * InvGamma:
			if s > 25:
				summe += sval[1] * (np.log(math.exp(math.log(w * t) + (-1 - s) * math.log(1 + t)) + math.exp(
					math.log(1. - w) + math.log(2.) + ((s + a) / 2.) * math.log(b) + float(mp.log(
						mp.besselk(-s + a, 2 * math.sqrt(b)).real)) - sp.gammaln(s + 1) - sp.gammaln(a))))
			else:
				try:
					summe += sval[1] * (np.log(math.exp(math.log(w * t) + (-1 - s) * math.log(1 + t)) + math.exp(
						math.log(1. - w) + math.log(2.) + ((s + a) / 2.) * math.log(b) + math.log(
							sp.kv(-s + a, 2 * math.sqrt(b))) - sp.gammaln(s + 1) - sp.gammaln(a))))
				except:
					summe += sval[1] * (np.log(math.exp(math.log(w * t) + (-1 - s)*math.log(1 + t)) + math.exp(math.log(1.-w) + math.log(2.) + (
						(s + a)/2.)*math.log(b) + float(mp.log(mp.besselk(-s + a, 2*math.sqrt(b)).real)) - sp.gammaln(s + 1) - sp.gammaln(a))))

	elif modC == 5:
		for sval in genes_by_sobs:
			s = sval[0]
			# *************** lambda ~ w * Gamma + (1-w) * Gamma:
			summe += sval[1] * (np.log(np.exp(
				np.log(w) + s * np.log(b) + (-s - a) * np.log(1 + b) + sp.gammaln(s + a) - sp.gammaln(
					s + 1) - sp.gammaln(a)) + np.exp(
				np.log(1. - w) + s * np.log(d) + (-s - g) * np.log(1 + d) + sp.gammaln(s + g) - sp.gammaln(
					s + 1) - sp.gammaln(g))))
	elif modC == 6:
		for sval in genes_by_sobs:
			s = sval[0]
			# *************** lambda ~ w * Gamma + (1-w) * InvGamma:
			if s > 25:
				summe += sval[1] * (np.log((w * math.exp(
					s * math.log(b) + (-s - a) * math.log(1 + b) + sp.gammaln(s + a) - sp.gammaln(s + 1) - sp.gammaln(
						a))) + ((1. - w) * math.exp(math.log(2.) + ((s + g) / 2.) * math.log(d) + float(mp.log(
					mp.besselk(-s + g, 2 * math.sqrt(d)).real)) - sp.gammaln(s + 1) - sp.gammaln(g)))))
			else:
				try:
					summe += sval[1] * (np.log((w * b ** s * (1 + b) ** (-s - a) * math.exp(
						sp.gammaln(s + a) - sp.gammaln(s + 1) - sp.gammaln(a))) + ((1. - w) * math.exp(
						math.log(2.) + ((s + g) / 2.) * math.log(d) + math.log(
							sp.kv(-s + g, 2 * math.sqrt(d))) - sp.gammaln(s + 1) - sp.gammaln(g)))))
				except:
					summe += sval[1] * (np.log((w * math.exp(s*math.log(b) + (-s-a)*math.log(1 + b) + sp.gammaln(s + a) - sp.gammaln(s+1) - sp.gammaln(a))) + (
						(1.-w) * math.exp(math.log(2.) + ((s + g)/2.)*math.log(d) + float(mp.log(mp.besselk(-s + g, 2*math.sqrt(d)).real)) - sp.gammaln(s + 1) - sp.gammaln(g)))))

	lastres = -summe
	if lastres > 1e8:
		return 0
	return -summe


sys.stderr.write("Use input columns:\n")
sys.stderr.write("gene_name\tmutation_effect\tmutated_base\tcontext\n")

# ************************************** COMMAND LINE ARGS ***************************************

infile = str(sys.argv[1])  # somatic mutation data input file
filepath = str(sys.argv[2])  # path to auxiliary input files folder
c_mode = int(sys.argv[3])  # 0 = trinucleotides, 1 = pentanucleotides
# mod_C		= int(sys.argv[4])		#	model choice: 0=all, 1=G, 2=IG, 3=EmixG, 4=EmixIG, 5=GmixG, 6=GmixIG
outname = str(sys.argv[4])  # name for the output file containing the q values

# ***************************** GLOBAL DEFINITIONS & AUXILIARY DATA ******************************

# muttypes = ["missense", "nonsense", "coding-synon", "intron", "utr-3", "utr-5", "IGR"]
triplets = ["AAA", "AAC", "AAG", "AAT", "CAA", "CAC", "CAG", "CAT", "GAA", "GAC", "GAG", "GAT", "TAA", "TAC", "TAG",
			"TAT", "ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC",
			"TCG", "TCT", "AGA", "AGC", "AGG", "AGT", "CGA", "CGC", "CGG", "CGT", "GGA", "GGC", "GGG", "GGT", "TGA",
			"TGC", "TGG", "TGT", "ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT",
			"TTA", "TTC", "TTG", "TTT"]
triplets_user = ["AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC",
				 "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT",
				 "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC",
				 "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT",
				 "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"]
bases = ["A", "C", "G", "T", "N"]
mod_choice_short = ["", "G", "IG", "EmixG", "EmixIG", "GmixG", "GmixIG"]  # model choice
mod_choice = ["", "Gamma(a,b): [a,b] =", "InverseGamma(a,b): [a,b] =", "w Exp(t) + (1-w) Gamma(a,b): [a,b,t,w] =",
			  "w Exp(t) + (1-w) InverseGamma(a,b): [a,b,t,w] =", "w Gamma(a,b) + (1-w) Gamma(g,d): [a,b,g,d,w] =",
			  "w Gamma(a,b) + (1-w) InverseGamma(g,d): [a,b,g,d,w] ="]  # model choice

# rep_no might be adjusted for the sake of a fastest run; just in case it becomes too heavy
rep_no = 50  # No. of independent runs to estimate model parameters (maximizing log-likelihood)


quintuplets = import_quintuplets("%s/quintuplets.txt" % filepath)
cancer_genes = import_special_genes("%s/COSMIC_genes_v80.txt" % filepath)
essential_genes = import_special_genes("%s/Wang_cell_essential_genes.txt" % filepath)
zero_genes = import_special_genes("%s/zero_genes.txt" % filepath)
if c_mode == 0:
	nuc_context_occs = import_context_freqs("%s/trinucleotide_occurrences_exons.txt" % filepath)
else:
	nuc_context_occs = import_context_freqs("%s/pentanucleotide_occurrences_exons.txt" % filepath)

# ************************************************************************************************************
# ************************************************************************************************************

sys.stderr.write("Running data preparation.\n")
sys.stderr.write("Use input columns:\n")
sys.stderr.write("gene_name\tmutation_effect\tmutated_base\tcontext\n")

#	(1)	Import mutation annotation file (maf) including header line, "context" is 0-based.
#	Format: ["gene", "muttype", "mutbase", "context"]
mutations = import_maf_data(infile, c_mode)
lengths = [[k, len(list(g))] for k, g in
		   it.groupby(sorted(mutations, key=lambda arg: arg["muttype"]), key=lambda arg: arg["muttype"])]
sys.stderr.write("%i SNVs imported.\n" % (len(mutations)))
for el in lengths:
	sys.stderr.write("%s\t%i\n" % (el[0], el[1]))

#	(1.1) Select the model based on the number of mutations
#	Specification supplied by Donate Weghorn
#	model choice: 1=G, 2=IG, 3=EmixG, 4=EmixIG, 5=GmixG, 6=GmixIG
# mod_C = 0
l = len(mutations)
if l < 12000:
	mod_C = 2
elif 12000 <= l <= 65000:
	mod_C = 4
else:
	mod_C = 6

repeat = True
while repeat:
	try:
		#	(2)	Import trinucleotides and codons for each gene.
		codons_by_gene = import_codons_by_gene("%s/codons_by_gene.txt.gz" % filepath)
		sys.stderr.write("Derive expected counts for %i genes.\n" % len(codons_by_gene))

		#	(3)	Build array with neutral background mutations for construction of neutral mutation matrix.
		neutral_mutations = [mut for mut in mutations if ["missense", "nonsense", "coding-synon", "utr-3", "utr-5"].count(
			mut["muttype"]) and cancer_genes.count(mut["gene"]) == 0 and essential_genes.count(mut["gene"]) == 0]
		sys.stderr.write(
			"Total number of coding mutations used for generation of mutation matrix, exluding likely selected genes: %i.\n" % len(
				neutral_mutations))
		neutral_muts_by_context = [[mut["context"], bases.index(
			mut["mutbase"])] for mut in neutral_mutations]

		#	(4)	Derive expected and observed mutation counts for all three categories per gene.
		res = export_expected_observed_mks_per_gene(codons_by_gene, mutations, neutral_muts_by_context, nuc_context_occs,
													c_mode)

		sys.stderr.write("Finished data preparation.\n")

		# ************************************************************************************************************

		sys.stderr.write("Running parameter estimation.\n")

		mks_type = []
		for gene in res:
			mks_type.append({"gene": gene[0], "exp": [float(gene[1]), float(gene[2]), float(gene[3])],
							"obs": [float(gene[4]), float(gene[5]), float(gene[6])], "len": int(gene[7])})
		mks_type = sorted(mks_type, key=lambda arg: arg["gene"])

		sys.stderr.write("Running ML routine %i times.\n" % rep_no)
		if mod_C > 6:
			sys.stderr.write("Not a valid model choice.\n")
			sys.exit()
		elif mod_C == 0:
			sys.stderr.write("Testing all six models.\n")
		else:
			sys.stderr.write("lam_s ~ %s.\n" % mod_choice_short[mod_C])
		sys.stderr.write("%i genes in total.\n" % len(mks_type))
		mks_type = [mks for mks in mks_type if (len(mks["gene"]) > 2 and mks["gene"][:2] == "OR" and [
													'0', '1', '2', '3', '4', '5', '6', '7', '8',
																		'9'].count(mks["gene"][2])) == 0]
		sys.stderr.write("Filtered out OR genes, leaving %i genes.\n" % len(mks_type))
		mks_type = [gene for gene in mks_type if zero_genes.count(gene["gene"]) == 0]
		sys.stderr.write(
			"Filtered out likely undercalled genes, leaving %i genes.\n" % len(mks_type))

		fout = open("Output/output_data_preparation_%s.txt" % outname, "w")
		# Output format: [gene, lm, lk, ls, mobs, kobs, sobs, Lgene]
		for gene in mks_type:
			fout.write("%s\t%f\t%f\t%f\t%i\t%i\t%i\t%i\n" % (
			gene["gene"], gene["exp"][0], gene["exp"][1], gene["exp"][2], gene["obs"][0], gene["obs"][1], gene["obs"][2],
			gene["len"]))
		fout.close()

		if mod_C == 1 or mod_C == 0:
			sys.stderr.write("Fitting model 1...\n")
			low_b = [1e-5 * random.uniform(1., 3.) for i in range(2)]
			up_b = [50. * random.uniform(1., 2.) for i in range(2)]
			cur_min_res = [0, 0, 1e20]
			for rep in range(rep_no):
				sys.stderr.write("%.f%% done.\r" % (100.*rep/rep_no))
				p_res = minimize_neg_ln_L([random.uniform(0.02, 10.), random.uniform(0.02, 10.)], neg_ln_L, mks_type, 1,
										[(low_b[0], up_b[0]), (low_b[1], up_b[1])], 2)
				if p_res[2] > 0 and p_res[2] < cur_min_res[2]:
					cur_min_res = p_res[:]
			if cur_min_res[2] == 1e20:
				sys.stderr.write("Could not find a converging solution for model 1.\n")
				if mod_C != 0:
					raise Exception
			fout = open("Output/param_estimates_%s_1.txt" % outname, "w")
			fout.write("%e, %e, %f, %i\n" % (cur_min_res[0], cur_min_res[1], cur_min_res[2], 1))
			fout.close()

		if mod_C == 2 or mod_C == 0:
			sys.stderr.write("Fitting model 2...\n")
			low_b = [1e-5 * random.uniform(1., 3.) for i in range(2)]
			up_b = [50. * random.uniform(1., 2.) for i in range(2)]
			cur_min_res = [0, 0, 1e20]
			for rep in range(rep_no):
				sys.stderr.write("%.f%% done.\r" % (100.*rep/rep_no))
				p_res = minimize_neg_ln_L([random.uniform(0.02, 10.), random.uniform(0.02, 10.)], neg_ln_L, mks_type, 2,
										[(low_b[0], up_b[0]), (low_b[1], up_b[1])], 2)
				if p_res[2] > 0 and p_res[2] < cur_min_res[2]:
					cur_min_res = p_res[:]
			if cur_min_res[2] == 1e20:
				sys.stderr.write("Could not find a converging solution for model 2.\n")
				if mod_C != 0:
					raise Exception
			fout = open("Output/param_estimates_%s_2.txt" % outname, "w")
			fout.write("%e, %e, %f, %i\n" % (cur_min_res[0], cur_min_res[1], cur_min_res[2], 2))
			fout.close()

		if mod_C == 3 or mod_C == 0:
			sys.stderr.write("Fitting model 3...\n")
			low_b = [1e-5 * random.uniform(1., 3.) for i in range(4)]
			up_b = [50. * random.uniform(1., 2.) for i in range(4)]
			cur_min_res = [0, 0, 0, 0, 1e20]
			for rep in range(rep_no):
				sys.stderr.write("%.f%% done.\r" % (100.*rep/rep_no))
				p_res = minimize_neg_ln_L([random.uniform(0.02, 10.), random.uniform(0.02, 10.), random.uniform(0.02, 10.),
										random.uniform(2e-5, 0.95)], neg_ln_L, mks_type, 3,
										[(low_b[0], up_b[0]), (low_b[1], up_b[1]), (low_b[2], up_b[2]), (low_b[3], 0.9999)],
										4)
				if p_res[4] > 0 and p_res[4] < cur_min_res[4]:
					cur_min_res = p_res[:]
			if cur_min_res[4] == 1e20:
				sys.stderr.write("Could not find a converging solution for model 3.\n")
				if mod_C != 0:
					raise Exception
			fout = open("Output/param_estimates_%s_3.txt" % outname, "w")
			fout.write("%e, %e, %e, %e, %f, %i\n" % (
			cur_min_res[0], cur_min_res[1], cur_min_res[2], cur_min_res[3], cur_min_res[4], 3))
			fout.close()

		if mod_C == 4 or mod_C == 0:
			sys.stderr.write("Fitting model 4...\n")
			low_b = [1e-5 * random.uniform(1., 3.) for i in range(4)]
			up_b = [50. * random.uniform(1., 2.) for i in range(4)]
			cur_min_res = [0, 0, 0, 0, 1e20]
			for rep in range(rep_no):
				sys.stderr.write("%.f%% done.\r" % (100.*rep/rep_no))
				p_res = minimize_neg_ln_L([random.uniform(0.02, 10.), random.uniform(0.02, 10.), random.uniform(0.02, 10.),
										random.uniform(2e-5, 0.95)], neg_ln_L, mks_type, 4,
										[(low_b[0], up_b[0]), (low_b[1], up_b[1]), (low_b[2], up_b[2]), (low_b[3], 0.9999)],
										4)
				if p_res[4] > 0 and p_res[4] < cur_min_res[4]:
					cur_min_res = p_res[:]
			if cur_min_res[4] == 1e20:
				sys.stderr.write("Could not find a converging solution for model 4.\n")
				if mod_C != 0:
					raise Exception
			fout = open("Output/param_estimates_%s_4.txt" % outname, "w")
			fout.write("%e, %e, %e, %e, %f, %i\n" % (
				cur_min_res[0], cur_min_res[1], cur_min_res[2], cur_min_res[3], cur_min_res[4], 4))
			fout.close()

		if mod_C == 5 or mod_C == 0:
			sys.stderr.write("Fitting model 5...\n")
			low_b = [1e-5 * random.uniform(1., 3.) for i in range(5)]
			up_b = [50. * random.uniform(1., 2.) for i in range(5)]
			cur_min_res = [0, 0, 0, 0, 0, 1e20]
			for rep in range(int(rep_no)):
				sys.stderr.write("%.f%% done.\r" % (100.*rep/rep_no))
				p_res = minimize_neg_ln_L(
					[random.uniform(0.02, 10.), random.uniform(0.02, 5.), random.uniform(0.02, 10.), random.uniform(0.02, 10.),
					random.uniform(2e-5, 0.95)], neg_ln_L, mks_type, 5, 
					[(low_b[0], up_b[0]), (low_b[1], up_b[1]), (low_b[2], up_b[2]), (low_b[3], up_b[3]), (low_b[4], 0.9999)], 5)
				if p_res[5] > 0 and p_res[5] < cur_min_res[5]:
					cur_min_res = p_res[:]
			if cur_min_res[5] == 1e20:
				sys.stderr.write("Could not find a converging solution for model 5.\n")
				if mod_C != 0:
					raise Exception
			fout = open("Output/param_estimates_%s_5.txt" % outname, "w")
			fout.write("%e, %e, %e, %e, %e, %f, %i\n" % (
			cur_min_res[0], cur_min_res[1], cur_min_res[2], cur_min_res[3], cur_min_res[4], cur_min_res[5], 5))
			fout.close()

		if mod_C == 6 or mod_C == 0:
			sys.stderr.write("Fitting model 6...\n")
			low_b = [1e-5 * random.uniform(1., 3.) for i in range(5)]
			up_b = [50. * random.uniform(1., 2.) for i in range(5)]
			cur_min_res = [0, 0, 0, 0, 0, 1e20]
			for rep in range(int(2 * rep_no)):
				sys.stderr.write("%.f%% done.\r" % (100.*rep/(2*rep_no)))
				p_res = minimize_neg_ln_L(
					[random.uniform(0.02, 10.), random.uniform(0.02, 5.), random.uniform(0.02, 10.), random.uniform(0.02, 10.),
					random.uniform(2e-5, 0.95)], neg_ln_L, mks_type, 6, 
					[(low_b[0], up_b[0]), (low_b[1], up_b[1]), (low_b[2], up_b[2]), (low_b[3], up_b[3]), (low_b[4], 0.9999)], 5)
				if p_res[5] > 0 and p_res[5] < cur_min_res[5]:
					cur_min_res = p_res[:]
			if cur_min_res[5] == 1e20:
				sys.stderr.write("Could not find a converging solution for model 6.\n")
				if mod_C != 0:
					raise Exception
			fout = open("Output/param_estimates_%s_6.txt" % outname, "w")
			fout.write("%e, %e, %e, %e, %e, %f, %i\n" % (
			cur_min_res[0], cur_min_res[1], cur_min_res[2], cur_min_res[3], cur_min_res[4], cur_min_res[5], 6))
			fout.close()
		
		repeat = False
	except Exception:
		if mod_C != 1:
			sys.stderr.write("Error using model %i trying with model 1\n" % mod_C)
			mod_C = 1 
			repeat = True           
		else:
			repeat = False

#!/usr/bin/env python

# ************************************************************
# * Cancer Bayesian Selection Estimation (CBaSE):			 *
# * Code accompanying Weghorn & Sunyaev, Nat. Genet. (2017). *
# *															 *
# * Author:		Donate Weghorn								 *
# *															 *
# * Copyright:	(c) 2017 Donate Weghorn						 *
# *															 *
# * License:		Public Domain							 *
# *															 *
# * Version:		Modified version built on Version 1.0 	 *
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
import gzip
import numpy as np
import random
import mpmath as mp
import scipy.stats as st
import itertools as it


# ************************************************************************************************************
#	FUNCTION DEFINITIONS

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
        return [[1, 1, 1, 1], [1, 0, 0, 0], [1, 0, 0, 0]]
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
    with gzip.open(filename) as fin:
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
    for pat in neut_cont_array:
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
                            cur_pent = codons_gene[cod_cnt - 1][2] + codon[2][:2]
                    elif i == 1:
                        cur_pent = codon[0] + codon[2][1:]
                    else:
                        if cod_cnt == len(codons_gene) - 1:
                            continue
                        else:
                            cur_pent = codon[1] + codons_gene[cod_cnt + 1][0][1:]
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

    summe = 0.
    if modC == 1:
        for gind in range(len(genes)):
            s = int(genes[gind]["obs"][2])
            # *************** lambda ~ Gamma:
            summe += (s * np.log(b) + (-s - a) * np.log(1 + b) + sp.gammaln(s + a) - sp.gammaln(s + 1) - sp.gammaln(a))
    elif modC == 2:
        for gind in range(len(genes)):
            s = int(genes[gind]["obs"][2])
            # *************** lambda ~ IG:
            if s > 25:
                summe += math.log(2.) + ((s + a) / 2.) * math.log(b) + mp.log(
                    mp.besselk(-s + a, 2 * math.sqrt(b))) - sp.gammaln(s + 1) - sp.gammaln(a)
            else:
                summe += math.log(2.) + ((s + a) / 2.) * math.log(b) + math.log(
                    sp.kv(-s + a, 2 * math.sqrt(b))) - sp.gammaln(s + 1) - sp.gammaln(a)
    elif modC == 3:
        for gind in range(len(genes)):
            s = int(genes[gind]["obs"][2])
            # *************** lambda ~ w * Exp + (1-w) * Gamma:
            summe += np.log(math.exp(math.log(w * t) + (-1 - s) * math.log(1 + t)) + math.exp(
                math.log(1. - w) + s * math.log(b) + (-s - a) * math.log(1 + b) + sp.gammaln(s + a) - sp.gammaln(
                    s + 1) - sp.gammaln(a)))
    elif modC == 4:
        for gind in range(len(genes)):
            s = int(genes[gind]["obs"][2])
            # *************** lambda ~ w * Exp + (1-w) * InvGamma:
            if s > 25:
                summe += np.log(math.exp(math.log(w * t) + (-1 - s) * math.log(1 + t)) + math.exp(
                    math.log(1. - w) + math.log(2.) + ((s + a) / 2.) * math.log(b) + mp.log(
                        mp.besselk(-s + a, 2 * math.sqrt(b))) - sp.gammaln(s + 1) - sp.gammaln(a)))
            else:
                summe += np.log(math.exp(math.log(w * t) + (-1 - s) * math.log(1 + t)) + math.exp(
                    math.log(1. - w) + math.log(2.) + ((s + a) / 2.) * math.log(b) + math.log(
                        sp.kv(-s + a, 2 * math.sqrt(b))) - sp.gammaln(s + 1) - sp.gammaln(a)))
    elif modC == 5:
        for gind in range(len(genes)):
            s = int(genes[gind]["obs"][2])
            # *************** lambda ~ w * Gamma + (1-w) * Gamma:
            summe += np.log(np.exp(
                np.log(w) + s * np.log(b) + (-s - a) * np.log(1 + b) + sp.gammaln(s + a) - sp.gammaln(
                    s + 1) - sp.gammaln(a)) + np.exp(
                np.log(1. - w) + s * np.log(d) + (-s - g) * np.log(1 + d) + sp.gammaln(s + g) - sp.gammaln(
                    s + 1) - sp.gammaln(g)))
    elif modC == 6:
        for gind in range(len(genes)):
            s = int(genes[gind]["obs"][2])
            # *************** lambda ~ w * Gamma + (1-w) * InvGamma:
            if s > 25:
                summe += np.log((w * math.exp(
                    s * math.log(b) + (-s - a) * math.log(1 + b) + sp.gammaln(s + a) - sp.gammaln(s + 1) - sp.gammaln(
                        a))) + ((1. - w) * math.exp(math.log(2.) + ((s + g) / 2.) * math.log(d) + mp.log(
                    mp.besselk(-s + g, 2 * math.sqrt(d))) - sp.gammaln(s + 1) - sp.gammaln(g))))
            else:
                summe += np.log((w * b ** s * (1 + b) ** (-s - a) * math.exp(
                    sp.gammaln(s + a) - sp.gammaln(s + 1) - sp.gammaln(a))) + ((1. - w) * math.exp(
                    math.log(2.) + ((s + g) / 2.) * math.log(d) + math.log(
                        sp.kv(-s + g, 2 * math.sqrt(d))) - sp.gammaln(s + 1) - sp.gammaln(g))))

    lastres = -summe
    if lastres > 1e8:
        sys.stderr.write(
            "-ln(L) too large. Run again or try changing s threshold in neg_ln_L() or omitting genes with large s.\n")
        return 0
    return -summe


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
            return (L * b) ** s * (1 + L * b) ** (-s - a) * math.gamma(s + a) / (math.gamma(s + 1) * math.gamma(a))

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
    for gene in genes:

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
            meant2 = int(ratm * sobs) * 2  # = E[m] * 2
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
                print "m_ppos outside limits:"
                print large_flag, m_ppos
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
                print "k_ppos outside limits:"
                print large_flag, k_ppos
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
                m_ppos = m_ppos.real
                m_pneg = m_pneg.real
                if m_ppos < 0 or m_ppos > 1:
                    cum_p = 0.
                    for x in range(mobs):
                        cum_p += pofx_given_s(x, sobs, L, ratm, 1)
                    m_pneg = cum_p + pofx_given_s(mobs, sobs, L, ratm, 1)
                    m_ppos = 1. - cum_p
                    m_ppos = m_ppos.real
                    m_pneg = m_pneg.real
                    if m_ppos < 0 or m_ppos > 1 or math.isnan(m_ppos):
                        sys.stderr.write("Setting p_m^pos --> 0 on gene %s (was %e).\n" % (gene["gene"], m_ppos))
                        m_ppos = 0.
                cum_p = 0.
                for x in range(kobs):
                    cum_p += pofx_given_s(x, sobs, L, ratk, large_flag)
                k_pneg = cum_p + pofx_given_s(kobs, sobs, L, ratk, large_flag)
                k_ppos = 1. - cum_p
                k_ppos = k_ppos.real
                k_pneg = k_pneg.real
                if k_ppos < 0 or k_ppos > 1 or math.isnan(k_ppos):
                    cum_p = 0.
                    for x in range(kobs):
                        cum_p += pofx_given_s(x, sobs, L, ratk, 1)
                    k_pneg = cum_p + pofx_given_s(kobs, sobs, L, ratk, 1)
                    k_ppos = 1. - cum_p
                    k_ppos = k_ppos.real
                    k_pneg = k_pneg.real
                    if k_ppos < 0 or k_ppos > 1:
                        sys.stderr.write("Setting p_k^pos --> 0 on gene %s (was %e).\n" % (gene["gene"], k_ppos))
                        k_ppos = 0.

            pvals.append([gene["gene"], m_pneg, k_pneg, m_ppos, k_ppos, pofx_given_s(0, sobs, L, ratm, 0),
                          pofx_given_s(0, sobs, L, ratk, 0), mobs, kobs, sobs])

    return pvals


def construct_histogram(var_array, bin_var):
    var_max = max(var_array) + bin_var
    hist = [0. for i in range(int(var_max / bin_var))]
    for var in var_array:
        hist[int(var / bin_var)] += 1
    return [[(i + 0.5) * bin_var, hist[i] / len(var_array)] for i in range(len(hist))]


def compute_phi_sim(pvals_array, ind1, ind2):
    all_phi = []
    for gene in pvals_array:
        if gene[ind1] == 0. or gene[ind2] == 0.:
            all_phi.append(float("inf"))
        else:
            all_phi.append(-math.log(gene[ind1]) - math.log(gene[ind2]))
    return all_phi


def compute_phi_obs(pvals_array, ind1, ind2):
    all_phi = []
    for gene in pvals_array:
        if gene[ind1] == 0. or gene[ind2] == 0.:
            cur_phi = float("inf")
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

    phi_sim_max = max(phi_sim_array)

    phi_pvals_obs = []
    for gene in sorted_phi:
        if gene["phi"] > phi_sim_max + 0.5 * bin_phi:
            phi_pvals_obs.append({"gene": gene["gene"], "p_phi": 0, "phi": gene["phi"], "mks": gene["mks"]})
        else:
            if abs(gene["phi"]) < 1e-30:
                phi_pvals_obs.append(
                    {"gene": gene["gene"], "p_phi": 1. - cumprob, "phi": gene["p0m"] * gene["p0k"], "mks": gene["mks"]})
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

    phi_qvals_adj = []
    for g in range(len(phi_qvals_obs)):
        cur_min = min([el["q_phi"] for el in phi_qvals_obs[g:]])
        gene = phi_qvals_obs[g]
        phi_qvals_adj.append(
            {"gene": gene["gene"], "p_phi": gene["p_phi"], "q_adj": min(1., cur_min), "phi": gene["phi"],
             "mks": gene["mks"]})

    return phi_qvals_adj


# ************************************************************************************************************
#	COMMAND LINE ARGS

infile = str(sys.argv[1])  # somatic mutation data input file
filepath = str(sys.argv[2])  # path to auxiliary input files folder
c_mode = int(sys.argv[3])  # 0 = trinucleotides, 1 = pentanucleotides
# mod_C		= int(sys.argv[4])		#	model choice: 1=G, 2=IG, 3=EmixG, 4=EmixIG, 5=GmixG, 6=GmixIG

# ************************************************************************************************************
#	GLOBAL DEFINITIONS & AUXILIARY DATA

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
run_no = 100  # No. of simulation replicates used for computing FDR

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
l = len(mutations)
if l < 12000:
    mod_C = 2
elif 12000 <= l <= 65000:
    mod_C = 4
else:
    mod_C = 6

#	(2)	Import trinucleotides and codons for each gene.
codons_by_gene = import_codons_by_gene("%s/codons_by_gene.txt.gz" % filepath)
sys.stderr.write("Derive expected counts for %i genes.\n" % len(codons_by_gene))

#	(3)	Build array with neutral background mutations for construction of neutral mutation matrix.
neutral_mutations = [mut for mut in mutations if ["missense", "nonsense", "coding-synon", "utr-3", "utr-5"].count(
    mut["muttype"]) and cancer_genes.count(mut["gene"]) == 0 and essential_genes.count(mut["gene"]) == 0]
sys.stderr.write(
    "Total number of coding mutations used for generation of mutation matrix, exluding likely selected genes: %i.\n" % len(
        neutral_mutations))
neutral_muts_by_context = [[mut["context"], bases.index(mut["mutbase"])] for mut in neutral_mutations]

#	(4)	Derive expected and observed mutation counts for all three categories per gene.
res = export_expected_observed_mks_per_gene(codons_by_gene, mutations, neutral_muts_by_context, nuc_context_occs,
                                            c_mode)

sys.stderr.write("Finished data preparation.\n")

# ************************************************************************************************************

sys.stderr.write("Running parameter_estimation.\n")

mks_type = []
for gene in res:
    mks_type.append({"gene": gene[0], "exp": [float(gene[1]), float(gene[2]), float(gene[3])],
                     "obs": [float(gene[4]), float(gene[5]), float(gene[6])], "len": int(gene[7])})
mks_type = sorted(mks_type, key=lambda arg: arg["gene"])

sys.stderr.write("Running ML routine %i times.\n" % rep_no)
if mod_C == 0:
    sys.stderr.write("Model selection from all six.\n")
else:
    sys.stderr.write("lam_s ~ %s.\n" % mod_choice_short[mod_C])
sys.stderr.write("%i genes in total.\n" % len(mks_type))
mks_type = [mks for mks in mks_type if (
            len(mks["gene"]) > 2 and mks["gene"][:2] == "OR" and ['0', '1', '2', '3', '4', '5', '6', '7', '8',
                                                                  '9'].count(mks["gene"][2])) == 0]
sys.stderr.write("Filtered out OR genes, leaving %i genes.\n" % len(mks_type))
mks_type = [gene for gene in mks_type if zero_genes.count(gene["gene"]) == 0]
sys.stderr.write(
    "Filtered out genes with M+K+S=0 in TumorPortal pan-cancer analysis, leaving %i genes.\n" % len(mks_type))

fout = open("output_data_preparation.txt", "w")
# Output format: [gene, lm, lk, ls, mobs, kobs, sobs, Lgene]
for gene in mks_type:
    fout.write("%s\t%f\t%f\t%f\t%i\t%i\t%i\t%i\n" % (
    gene["gene"], gene["exp"][0], gene["exp"][1], gene["exp"][2], gene["obs"][0], gene["obs"][1], gene["obs"][2],
    gene["len"]))
fout.close()

if mod_C == 1:
    sys.stderr.write("Fitting model 1...\n")
    low_b = [1e-5 * random.uniform(1., 3.) for i in range(2)]
    up_b = [50. * random.uniform(1., 2.) for i in range(2)]
    cur_min_res = [0, 0, 1e20]
    for rep in range(rep_no):
        p_res = minimize_neg_ln_L([random.uniform(0.02, 10.), random.uniform(0.02, 10.)], neg_ln_L, mks_type, 1,
                                  [(low_b[0], up_b[0]), (low_b[1], up_b[1])], 2)
        if p_res[2] > 0 and p_res[2] < cur_min_res[2]:
            cur_min_res = p_res[:]
    if cur_min_res[2] == 1e20:
        sys.stderr.write("Could not find a converging solution.\n")
    fout = open("param_estimates_1.txt", "w")
    fout.write("%e, %e, %f, %i\n" % (cur_min_res[0], cur_min_res[1], cur_min_res[2], mod_C))
    fout.close()

elif mod_C == 2:
    sys.stderr.write("Fitting model 2...\n")
    low_b = [1e-5 * random.uniform(1., 3.) for i in range(2)]
    up_b = [50. * random.uniform(1., 2.) for i in range(2)]
    cur_min_res = [0, 0, 1e20]
    for rep in range(rep_no):
        p_res = minimize_neg_ln_L([random.uniform(0.02, 10.), random.uniform(0.02, 10.)], neg_ln_L, mks_type, 2,
                                  [(low_b[0], up_b[0]), (low_b[1], up_b[1])], 2)
        if p_res[2] > 0 and p_res[2] < cur_min_res[2]:
            cur_min_res = p_res[:]
    if cur_min_res[2] == 1e20:
        sys.stderr.write("Could not find a converging solution.\n")
    fout = open("param_estimates_2.txt", "w")
    fout.write("%e, %e, %f, %i\n" % (cur_min_res[0], cur_min_res[1], cur_min_res[2], mod_C))
    fout.close()

elif mod_C == 3:
    sys.stderr.write("Fitting model 3...\n")
    low_b = [1e-5 * random.uniform(1., 3.) for i in range(4)]
    up_b = [50. * random.uniform(1., 2.) for i in range(4)]
    cur_min_res = [0, 0, 0, 0, 1e20]
    for rep in range(rep_no):
        p_res = minimize_neg_ln_L([random.uniform(0.02, 10.), random.uniform(0.02, 10.), random.uniform(0.02, 10.),
                                   random.uniform(2e-5, 0.95)], neg_ln_L, mks_type, 3,
                                  [(low_b[0], up_b[0]), (low_b[1], up_b[1]), (low_b[2], up_b[2]), (low_b[3], 0.9999)],
                                  4)
        if p_res[4] > 0 and p_res[4] < cur_min_res[4]:
            cur_min_res = p_res[:]
    if cur_min_res[4] == 1e20:
        sys.stderr.write("Could not find a converging solution.\n")
    fout = open("param_estimates_3.txt", "w")
    fout.write("%e, %e, %e, %e, %f, %i\n" % (
    cur_min_res[0], cur_min_res[1], cur_min_res[2], cur_min_res[3], cur_min_res[4], mod_C))
    fout.close()

elif mod_C == 4:
    sys.stderr.write("Fitting model 4...\n")
    low_b = [1e-5 * random.uniform(1., 3.) for i in range(4)]
    up_b = [50. * random.uniform(1., 2.) for i in range(4)]
    cur_min_res = [0, 0, 0, 0, 1e20]
    for rep in range(rep_no):
        p_res = minimize_neg_ln_L([random.uniform(0.02, 10.), random.uniform(0.02, 10.), random.uniform(0.02, 10.),
                                   random.uniform(2e-5, 0.95)], neg_ln_L, mks_type, 4,
                                  [(low_b[0], up_b[0]), (low_b[1], up_b[1]), (low_b[2], up_b[2]), (low_b[3], 0.9999)],
                                  4)
        if p_res[4] > 0 and p_res[4] < cur_min_res[4]:
            cur_min_res = p_res[:]
    if cur_min_res[4] == 1e20:
        sys.stderr.write("Could not find a converging solution.\n")
    fout = open("param_estimates_4.txt", "w")
    fout.write("%e, %e, %e, %e, %f, %i\n" % (
    cur_min_res[0], cur_min_res[1], cur_min_res[2], cur_min_res[3], cur_min_res[4], mod_C))
    fout.close()

elif mod_C == 5:
    sys.stderr.write("Fitting model 5...\n")
    low_b = [1e-5 * random.uniform(1., 3.) for i in range(5)]
    up_b = [50. * random.uniform(1., 2.) for i in range(5)]
    cur_min_res = [0, 0, 0, 0, 0, 1e20]
    for rep in range(int(2 * rep_no)):
        p_res = minimize_neg_ln_L(
            [random.uniform(0.02, 10.), random.uniform(0.02, 5.), random.uniform(0.02, 10.), random.uniform(0.02, 10.),
             random.uniform(2e-5, 0.95)], neg_ln_L, mks_type, 5,
            [(low_b[0], up_b[0]), (low_b[1], up_b[1]), (low_b[2], up_b[2]), (low_b[3], up_b[3]), (low_b[4], 0.9999)], 5)
        if p_res[5] > 0 and p_res[5] < cur_min_res[5]:
            cur_min_res = p_res[:]
    if cur_min_res[5] == 1e20:
        sys.stderr.write("Could not find a converging solution.\n")
    fout = open("param_estimates_5.txt", "w")
    fout.write("%e, %e, %e, %e, %e, %f, %i\n" % (
    cur_min_res[0], cur_min_res[1], cur_min_res[2], cur_min_res[3], cur_min_res[4], cur_min_res[5], mod_C))
    fout.close()

elif mod_C == 6:
    sys.stderr.write("Fitting model 6...\n")
    low_b = [1e-5 * random.uniform(1., 3.) for i in range(5)]
    up_b = [50. * random.uniform(1., 2.) for i in range(5)]
    cur_min_res = [0, 0, 0, 0, 0, 1e20]
    for rep in range(int(2 * rep_no)):
        p_res = minimize_neg_ln_L(
            [random.uniform(0.02, 10.), random.uniform(0.02, 5.), random.uniform(0.02, 10.), random.uniform(0.02, 10.),
             random.uniform(2e-5, 0.95)], neg_ln_L, mks_type, 6,
            [(low_b[0], up_b[0]), (low_b[1], up_b[1]), (low_b[2], up_b[2]), (low_b[3], up_b[3]), (low_b[4], 0.9999)], 5)
        if p_res[5] > 0 and p_res[5] < cur_min_res[5]:
            cur_min_res = p_res[:]
    if cur_min_res[5] == 1e20:
        sys.stderr.write("Could not find a converging solution.\n")
    fout = open("param_estimates_6.txt", "w")
    fout.write("%e, %e, %e, %e, %e, %f, %i\n" % (
    cur_min_res[0], cur_min_res[1], cur_min_res[2], cur_min_res[3], cur_min_res[4], cur_min_res[5], mod_C))
    fout.close()

else:
    sys.stderr.write("Invalid choice for mod_C: %i.\n" % mod_C)
    sys.exit()

# ************************************************************************************************************

sys.stderr.write("Running p_values.\n")

params = cur_min_res[:-1]

sys.stderr.write("simulation_runs\t= %i\n" % (run_no))
sys.stderr.write("Computing simulated p-values.\n")
pvals_sim = compute_p_values(params, mks_type, [mod_C, 1, run_no])
sys.stderr.write("Computing real p-values.\n")
pvals_obs = compute_p_values(params, mks_type, [mod_C, 0, 1])
# [gene["gene"], m_pneg, k_pneg, m_ppos, k_ppos, pofx_given_s(0, sobs, L, ratm), pofx_given_s(0, sobs, L, ratk)]

#	Negative selection analysis
# phi_sim = compute_phi_sim(pvals_sim, 1, 2)
# gene_phi_obs = compute_phi_obs(pvals_obs, 1, 2)
# q_neg_adj = FDR_discrete(phi_sim, gene_phi_obs, 0.02, 0.000001)

#	Positive selection analysis
phi_sim = compute_phi_sim(pvals_sim, 3, 4)
gene_phi_obs = compute_phi_obs(pvals_obs, 3, 4)
q_pos_adj = FDR_discrete(phi_sim, gene_phi_obs, 0.02, 0.000001)

# q_neg_sorted = sorted(q_neg_adj, key=lambda arg: arg["gene"])
q_pos_sorted = sorted(q_pos_adj, key=lambda arg: arg["gene"])
pvals = [val["p_phi"] for val in q_pos_sorted]
qvals = multipletests(pvals, method='fdr_bh')[1]

fout = open("q_values_output.txt", "w")
fout.write("%s\t%s\n" % (mod_choice[mod_C], params))
# fout.write(
#   "gene\tp_phi_neg\tq_phi_neg\tphi_neg\tp_phi_pos\tq_phi_pos\tphi_pos_or_p(m=0|s)*p(k=0|s)\tm_obs\tk_obs\ts_obs\n")
fout.write(
    "gene\tp_pos\tq_pos\tmis_obs\tnon_obs\tsyn_obs\n")
for g in range(len(q_pos_sorted)):
    # fout.write("%s\t%e\t%e\t%f\t%e\t%e\t%f\t%i\t%i\t%i\n"
    #            % (q_neg_sorted[g]["gene"], q_neg_sorted[g]["p_phi"], q_neg_sorted[g]["q_adj"], q_neg_sorted[g]["phi"],
    #             q_pos_sorted[g]["p_phi"], q_pos_sorted[g]["q_adj"], q_pos_sorted[g]["phi"], q_neg_sorted[g]["mks"][0],
    #             q_neg_sorted[g]["mks"][1], q_neg_sorted[g]["mks"][2]))
    fout.write("%s\t%e\t%e\t%i\t%i\t%i\n"
               % (q_pos_sorted[g]["gene"], pvals[g], qvals[g], q_pos_sorted[g]["mks"][0],
                  q_pos_sorted[g]["mks"][1], q_pos_sorted[g]["mks"][2]))
fout.close()

## This script is to run the likelihood ratio test with Godambe adjustment for a specified taxa pair input as an argument.
## To run this for multiple taxa pairs, simply send in all pairs as arguments in a bash loop, running one pair per loop iteration.

import sys
import os
import csv
import dadi
import Models_2D
import dadi.Godambe

#get taxa pair as system argument
t = sys.argv[1]
m = sys.argv[2]

#import models
func_nomig = Models_2D.no_mig
func_ex_nomig = dadi.Numerics.make_extrap_log_func(func_nomig)

func_symmig = Models_2D.sym_mig
func_ex_symmig = dadi.Numerics.make_extrap_log_func(func_symmig)

func_ancsymmig = Models_2D.anc_sym_mig
func_ex_ancsymmig = dadi.Numerics.make_extrap_log_func(func_ancsymmig)

func_scsymmig = Models_2D.sec_contact_sym_mig
func_ex_scsymmig = dadi.Numerics.make_extrap_log_func(func_scsymmig)

func_CMG = Models_2D.CMG
func_ex_CMG = dadi.Numerics.make_extrap_log_func(func_CMG)

func_SCG = Models_2D.SCG
func_ex_SCG = dadi.Numerics.make_extrap_log_func(func_SCG)

#for t in taxa:
print(t)
print(m)

#import data SFS
if t == 'cppar_cpery':
    fs = dadi.Spectrum.from_file('../site_freq_spectra/sfs/2dsfs_cppar_cpery_p1_fold0_ancmol_genesum.sfs')
else:
    fs = dadi.Spectrum.from_file('../site_freq_spectra/sfs/2dsfs_%s_p1_fold0_genesum.sfs' % (t))

#fold SFS and get sample sizes
fs = fs.fold()
ns = fs.sample_sizes

#set grid points
pts = [sum(ns),sum(ns)+10,sum(ns)+20]

#import bootstraps
all_boot = [dadi.Spectrum.from_file('../site_freq_spectra/sfs/gl_boot_sfs/%s/boot%s.sfs' % (t,ii)) for ii in range(1,101)]
all_boot = [b.fold() for b in all_boot]

#read in optimized parameter values for no mig
with open('outputs/%s/gl_run/results_summary/%s_fold_no_mig_gl_results_summary.txt' % (t,t)) as f:
    reader = csv.reader(f, delimiter = '\t')
    no_mig = list(reader)
no_mig = no_mig[1:4]
no_mig.sort(key = lambda x: x[4])
no_mig = no_mig[0][7:10]
no_mig = [float(i) for i in no_mig]

#read in optimized parameter values for sym mig
with open('outputs/%s/gl_run/results_summary/%s_fold_sym_mig_gl_results_summary.txt' % (t,t)) as f:
    reader = csv.reader(f, delimiter = '\t')
    sym_mig = list(reader)
sym_mig = sym_mig[1:4]
sym_mig.sort(key = lambda x: x[4])
sym_mig = sym_mig[0][7:11]
sym_mig = [float(i) for i in sym_mig]

#read in optimized parameter values for anc mig
with open('outputs/%s/gl_run/results_summary/%s_fold_anc_sym_mig_gl_results_summary.txt' % (t,t)) as f:
    reader = csv.reader(f, delimiter = '\t')
    anc_sym_mig = list(reader)
anc_sym_mig = anc_sym_mig[1:4]
anc_sym_mig.sort(key = lambda x: x[4])
anc_sym_mig = anc_sym_mig[0][7:12]
anc_sym_mig = [float(i) for i in anc_sym_mig]

#read in optimized param values for sec contact
with open('outputs/%s/gl_run/results_summary/%s_fold_sec_contact_sym_mig_gl_results_summary.txt' % (t,t)) as f:
    reader = csv.reader(f, delimiter = '\t')
    sec_contact_sym_mig = list(reader)
sec_contact_sym_mig = sec_contact_sym_mig[1:4]
sec_contact_sym_mig.sort(key = lambda x: x[4])
sec_contact_sym_mig = sec_contact_sym_mig[0][7:12]
sec_contact_sym_mig = [float(i) for i in sec_contact_sym_mig]

#read in optimized parameter values for CMG
with open('outputs/%s/gl_run/results_summary/%s_fold_CMG_gl_results_summary.txt' % (t,t)) as f:
    reader = csv.reader(f, delimiter = '\t')
    CMG = list(reader)
CMG = CMG[1:4]
CMG.sort(key = lambda x: x[4])
CMG = CMG[0][7:13]
CMG = [float(i) for i in CMG]

#read in optimized parameter values for SCG
with open('outputs/%s/gl_run/results_summary/%s_fold_SCG_gl_results_summary.txt' % (t,t)) as f:
    reader = csv.reader(f, delimiter = '\t')
    SCG = list(reader)
SCG = SCG[1:4]
SCG.sort(key = lambda x: x[4])
SCG = SCG[0][7:14]
SCG = [float(i) for i in SCG]

#evaluate models with optimized parameters
ll_nomig = dadi.Inference.ll_multinom(func_ex_nomig(no_mig,ns,pts),fs)
ll_symmig = dadi.Inference.ll_multinom(func_ex_symmig(sym_mig,ns,pts),fs)
ll_ancmig = dadi.Inference.ll_multinom(func_ex_ancsymmig(anc_sym_mig,ns,pts),fs)
ll_scsymmig = dadi.Inference.ll_multinom(func_ex_scsymmig(sec_contact_sym_mig,ns,pts),fs)
ll_CMG = dadi.Inference.ll_multinom(func_ex_CMG(CMG,ns,pts),fs)
ll_SCG = dadi.Inference.ll_multinom(func_ex_SCG(SCG,ns,pts),fs)

#run adjusted LRT based on which model is input. simpler model is automatically selected.
if m == 'sym_mig':

    lrt_file = open('outputs/pop_pairs_lrt_sym_mig.tsv', 'a+')
    #lrt_file.write('pair\tll_nomig\tll_symmig\tpvalGIM\n')

    print('no mig vs. sym mig')

    print('likelihood of no mig: {}'.format(ll_nomig))
    print('likelihood of sym mig: {}'.format(ll_symmig))

    p_lrt_sym = list(no_mig)
    p_lrt_sym.insert(2,0)

    adj_sym = dadi.Godambe.LRT_adjust(func_ex_symmig, pts, all_boot, p_lrt_sym, fs, nested_indices=[2], multinom=True)
    D_sym = 2*(ll_symmig - ll_nomig)
    D_adj_sym = D_sym*adj_sym
    pval_adj_sym = dadi.Godambe.sum_chi2_ppf(D_adj_sym, weights=(0.5,0.5))

    print('p-value for rejecting no_mig with sym_mig (GIM adjust): {0:.4f}\n'.format(pval_adj_sym))
    lrt_file.write('%s\t%s\t%s\t%s\n' % (t,ll_nomig,ll_symmig,pval_adj_sym))
    lrt_file.close

elif m == 'anc_sym_mig':

    lrt_file = open('outputs/pop_pairs_lrt_anc_sym_mig.tsv', 'a+')
    #lrt_file.write('pair\tll_symmig\tll_ancsymmig\tpvalGIM\n')

    print('sym mig vs. anc mig')

    print('likelihood of symmetric migration model: {}'.format(ll_symmig))
    print('likelihood of ancestral migration model: {}'.format(ll_ancmig))

    p_lrt_anc = list(sym_mig)
    p_lrt_anc.insert(3,0)

    adj_anc = dadi.Godambe.LRT_adjust(func_ex_ancsymmig, pts, all_boot, p_lrt_anc, fs, nested_indices=[3], multinom=True)
    D_anc = 2*(ll_ancmig - ll_symmig)
    D_adj_anc = D_anc*adj_anc
    pval_adj_anc = dadi.Godambe.sum_chi2_ppf(D_adj_anc, weights=(0.5,0.5))

    print('p-value for rejecting sym_mig with anc_sym_mig (GIM adjust): {0:.4f}\n'.format(pval_adj_anc))
    lrt_file.write('%s\t%s\t%s\t%s\n' % (t,ll_symmig,ll_ancmig,pval_adj_anc))
    lrt_file.close

elif m == 'sec_contact_sym_mig':

    lrt_file = open('outputs/pop_pairs_lrt_sc_sym_mig.tsv', 'a+')
    #lrt_file.write('pair\tll_symmig\tll_scsymmig\tpvalGIM\n')

    print('sym mig vs. sec con')

    print('likelihood of symmetric migration model: {}'.format(ll_symmig))
    print('likelihood of secondary contact model: {}'.format(ll_scsymmig))

    p_lrt_sc = list(sym_mig)
    p_lrt_sc.insert(3,0)

    adj_sc = dadi.Godambe.LRT_adjust(func_ex_scsymmig, pts, all_boot, p_lrt_sc, fs, nested_indices=[3], multinom=True)
    D_sc = 2*(ll_scsymmig - ll_symmig)
    D_adj_sc = D_sc*adj_sc
    pval_adj_sc = dadi.Godambe.sum_chi2_ppf(D_adj_sc, weights=(0.5,0.5))

    print('p-value for rejecting sym_mig with sec_contact_sym_mig (GIM adjust): {0:.4f}\n'.format(pval_adj_sc))
    lrt_file.write('%s\t%s\t%s\t%s\n' % (t,ll_symmig,ll_scsymmig,pval_adj_sc))
    lrt_file.close

elif m == 'SCG':

    lrt_file = open('outputs/pop_pairs_lrt_SCG.tsv', 'a+')

    print('CMG vs. SCG')

    print('likelihood of CMG model: {}'.format(ll_CMG))
    print('likelihood of SCG model: {}'.format(ll_SCG))

    p_lrt_SCG = list(CMG)
    p_lrt_SCG.insert(5,0)

    print('CMG best fit params under SCG: {}'.format(p_lrt_SCG))
    print('likelihood of SCG model using best fit CMG params: {}'.format(dadi.Inference.ll_multinom(func_ex_SCG(p_lrt_SCG,ns,pts),fs)))

    adj_SCG = dadi.Godambe.LRT_adjust(func_ex_SCG, pts, all_boot, p_lrt_SCG, fs, nested_indices=[5], multinom=True)
    D_SCG = 2*(ll_SCG - ll_CMG)
    D_adj_SCG = D_SCG*adj_SCG
    pval_SCG = dadi.Godambe.sum_chi2_ppf(D_SCG, weights=(0.5,0.5))
    pval_adj_SCG = dadi.Godambe.sum_chi2_ppf(D_adj_SCG, weights=(0.5,0.5))

    print('p-value for rejecting CMG with SCG (GIM adjust): {0:.4f}\n'.format(pval_adj_SCG))
    lrt_file.write('%s\t%s\t%s\t%s\t%s\n' % (t,ll_CMG,ll_SCG,pval_SCG,pval_adj_SCG))
    lrt_file.close



else:
    print('not a valid model name')

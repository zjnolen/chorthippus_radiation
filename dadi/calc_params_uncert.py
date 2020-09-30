# As all population pair models converged on secondary contact, we evaluated parameter values under this model.

# We calculated uncertainty using the Godambe Information Matrix method implemented in dadi. This outputs parameter uncertainties for all parameters including theta. From these parameters, real values can be calculated and real uncertainties can be calculated by propogating uncertainty from the parameter uncertainties.

import sys
import os
import csv
import dadi
import Models_2D
import dadi.Godambe
import math

#Define taxa pairs, demographic model, and mutation rate
taxa = ['cppar_cpery','cbig_b_cbru_b','cbig_b_cbru_s','cbig_e_cbru_b','cbig_e_cbru_s','cmol_b_cbig_b','cmol_b_cbig_e','cmol_b_cbru_b','cmol_b_cbru_s','cmol_e_cbig_b','cmol_e_cbig_e','cmol_e_cbru_b','cmol_e_cbru_s','cbig_b_crub_a','cbig_e_crub_a','cbru_b_crub_a','cbru_s_crub_a','cmol_b_crub_a','cmol_e_crub_a']

func_scsymmig = Models_2D.sec_contact_sym_mig
func_ex_scsymmig = dadi.Numerics.make_extrap_log_func(func_scsymmig)

mu = 1

#Build file to output parameter values
param_file = open('outputs/pop_pairs_gl_params.tsv', 'w')
param_file.write('sfs_type\ttaxa\tnref\tn1\tn2\tm\ttsplit\ttsc\tnref_sd\tn1_sd\tn2_sd\tm_sd\ttsplit_sd\ttsc_sd\n')

#Loop process over all pairs
for t in taxa:
    #Read in SFS
    if t == 'cppar_cpery':
        fs = dadi.Spectrum.from_file('../site_freq_spectra/sfs/2dsfs_%s_p1_fold0_ancmol_genesum.sfs' % (t))
    else:
        fs = dadi.Spectrum.from_file('../site_freq_spectra/sfs/2dsfs_%s_p1_fold0_genesum.sfs' % (t))

    #Get sample size
    ns = fs.sample_sizes

    #Fold SFS
    fs = fs.fold()

    #Set grid points from sample size
    pts = [sum(ns),sum(ns)+10,sum(ns)+20]

    #Load in bootstraps
    all_boot = [dadi.Spectrum.from_file('../site_freq_spectra/sfs/gl_boot_sfs/%s/boot%s.sfs' % (t,ii)) for ii in range(1,101)]
    all_boot = [b.fold() for b in all_boot]

    #Read in optimized parameters
    with open('outputs/%s/gl_run/results_summary/%s_fold_sec_contact_sym_mig_gl_results_summary.txt' % (t,t)) as f:
        reader = csv.reader(f, delimiter = '\t')
        sec_contact_sym_mig = list(reader)
    sec_contact_sym_mig = sec_contact_sym_mig[1:4]
    sec_contact_sym_mig.sort(key = lambda x: x[4])
    sec_contact_sym_mig = sec_contact_sym_mig[0]
    sec_contact_sym_mig = map(sec_contact_sym_mig.__getitem__, [7,8,9,10,11,6])
    popt = [float(i) for i in sec_contact_sym_mig]

    print('Calculating parameter values and uncertainties for %s\n' % (t))

    print('Best fit model parameters under secondary contact model + theta: %s\n' % (popt))

    uncerts = dadi.Godambe.GIM_uncert(func_ex_scsymmig, pts, all_boot, popt[0:5], fs, multinom=True, return_GIM=True)

    print('Estimated parameter standard deviations from GIM: %s\n' % (uncerts[0]))

    # Get effective sequencing length
    if t == 'cppar_cpery':
        with open('../site_freq_spectra/sfs/2dsfs_%s_p1_fold0_ancmol_genesum.sfs' % (t)) as f:
            reader = csv.reader(f, delimiter = ' ')
            plain_sfs = list(reader)
    else:
        with open('../site_freq_spectra/sfs/2dsfs_%s_p1_fold0_genesum.sfs' % (t)) as f:
            reader = csv.reader(f, delimiter = ' ')
            plain_sfs = list(reader)
    plain_sfs = plain_sfs[1]
    plain_sfs = [float(i) for i in plain_sfs]
    L = sum(plain_sfs)
    print('Effective sequencing length: {0}\n'.format(L))

    # Calculate real params
    Nref = popt[5]/(4*mu*L)
    N1 = Nref * popt[0]
    N2 = Nref * popt[1]
    M = popt[2]/(2*Nref)
    T1 = (popt[3]*popt[5])/(2*mu*L)
    Tsc = (popt[4]*popt[5])/(2*mu*L)
    Tsplit = T1 + Tsc

    preal = [Nref,N1,N2,M,Tsplit,Tsc]

    print('Real parameter values (Nref, N1, N2, m, T1, T2): {0}\n'.format(preal))

    #Extrapolate uncertainty for each value

    import numpy.linalg
    covarmatrix = numpy.linalg.inv(uncerts[1])

    dtheta = (popt[3]+popt[4])/(2*L*mu)
    dt = popt[5]/(2*L*mu)

    Nref_sd = uncerts[0][5] * (1/(4*mu*L))
    N1_sd = N1*(math.sqrt((uncerts[0][5]/popt[5]) ** 2 + (uncerts[0][0]/popt[0]) ** 2 + (2*(covarmatrix[0][5]/(popt[5]*popt[0])))))
    N2_sd = N2*(math.sqrt((uncerts[0][5]/popt[5]) ** 2 + (uncerts[0][1]/popt[1]) ** 2 + (2*(covarmatrix[1][5]/(popt[5]*popt[1])))))
    M_sd = M*(math.sqrt((uncerts[0][5]/popt[5]) ** 2 + (uncerts[0][2]/popt[2]) ** 2 - (2*(covarmatrix[2][5]/(popt[5]*popt[2])))))
    T1_sd = T1*(math.sqrt((uncerts[0][5]/popt[5]) ** 2 + (uncerts[0][3]/popt[3]) ** 2 + (2*(covarmatrix[3][5]/(popt[5]*popt[3])))))
    Tsc_sd = Tsc*(math.sqrt((uncerts[0][5]/popt[5]) ** 2 + (uncerts[0][4]/popt[4]) ** 2 + (2*(covarmatrix[4][5]/(popt[5]*popt[4])))))
    Trealcov = ((1/(2*mu*L)) ** 2) * (popt[3]*popt[4]*(uncerts[0][5]**2) + (popt[5] ** 2) * covarmatrix[3][4] + popt[5] * popt[4] * covarmatrix[3][5] + popt[5] * popt[3] * covarmatrix[4][5])
    Tsplit_sd = math.sqrt((T1_sd ** 2) + (Tsc_sd ** 2) + (2 * Trealcov))

    print(Nref_sd,N1_sd,N2_sd,M_sd,Tsplit_sd,Tsc_sd)

    param_file.write('gl\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (t,Nref,N1,N2,M,Tsplit,Tsc,Nref_sd,N1_sd,N2_sd,M_sd,Tsplit_sd,Tsc_sd))

param_file.close()

import sys
import os
from datetime import datetime
import time
import matplotlib
matplotlib.use('agg')
import numpy
import dadi
import pylab
import Optimize_Functions
import Models_2D

#In this script, a single demographic model will be tested for one pair of taxa. This performs one replicate of a full optimization run. It will add the optimized parameter and maximum likelihood to compiled output file for this taxa pair and model. This script was written with

#argument order: taxa_pair demo_model run_number fold/unfold sfs_type

#taxa_pair - The two populations to be analyzed, in the format of their sfs file name
#demo_model - The demographic model to fit the data to, leave out the 'Models_2D.'
#run_number - Number to append to the output names to differentiate parallel runs.
#fold/unfold - Whether the SFS should be folded or not. Use 'fold' or 'unfold'
#sfs_type - Whether to use the spectra produced by the full genotype likelihood data, or use the first bootstrap which is sampled at one site per gene. Can only take values 'gl' or 'snp'. Ultimately, we did not use the snp option in this manuscript, so it has been depracated and commented out.

taxa = sys.argv[1]
model = sys.argv[2]
run_num = int(sys.argv[3])
is_folded = sys.argv[4]
#sfs_type = sys.argv[5]
t1 = taxa[:(len(taxa)/2)]
t2 = taxa[-(len(taxa)/2):]

#Set output directory and make it if it doesn't exist.

if not run_num == 1:
	time.sleep(5)

out_dir = "outputs/%s/%s_run" % (taxa,sfs_type)

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

if not os.path.exists("%s/results_summary" % (out_dir)):
    os.makedirs("%s/results_summary" % (out_dir))

#Set temp log file directory and make it if it doesn't exist.
temp_dir = "temp/%s" % (out_dir)

if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)

#Import sfs. Will import based on whether gl or snp SFS was selected in the input arguments.
print "\n\nImporting frequency spectra..."

if taxa == 'cppar_cpery':
		fs_gl = dadi.Spectrum.from_file("../site_freq_spectra/sfs/2dsfs_%s_p1_fold0_ancmol_genesum.sfs" % (taxa))
else:
		fs_gl = dadi.Spectrum.from_file("../site_freq_spectra/sfs/2dsfs_%s_p1_fold0_genesum.sfs" % (taxa))

ns = fs_gl.sample_sizes

fs = fs_gl

#If returning to using snp SFS, uncomment lines below and comment out line above.
#if sfs_type == 'gl':
#	fs = fs_gl
#elif sfs_type == 'snp':
#	fs = dadi.Misc.make_data_dict("../genotype_calling/outputs/bootstraps/%s/boot_0.snp" % (taxa))
#	fs = dadi.Spectrum.from_data_dict(fs, [t1,t2], ns, polarized=True)
#else:
#	print "Please use 'gl' or 'snp' to describe the data type you would like to fit the model to. See script for more information on input arguments."


#Show SFS to double check
print "Site Freq Spectrum:"

print fs

#Fold if required
if is_folded == "fold":
	print "\n\nFolding frequency spectra..."
	fs = fs.fold()
elif is_folded == "unfold":
	print "\n\nContinuing with unfolded spectra..."
else:
	print "\n\nSyntax incorrect for folding argument, please use 'fold' or 'unfold'"

#Set priors that change based on spectrum type.
#if sfs_type == 'gl':

	initials = {
		"no_mig":[2,2,0.5],
		"sym_mig":[2,2,0.5,0.5],
		"CMG":[1,1,5,5,0.5,1],
		"anc_sym_mig":[2,2,1,0.5,0.5],
		"sec_contact_sym_mig":[2,2,1,0.5,0.5],
		"SCG":[1,1,5,5,0.5,1,1],
		}

	uppers = {
		"no_mig":[6,6,2],
		"sym_mig":[5,5,2,2],
		"CMG":[None,None,None,None,None,None],
		"anc_sym_mig":[6,6,2,2,2],
		"sec_contact_sym_mig":[5,5,5,2,2],
		"SCG":[100,100,100,100,5,10,10],
		}

# elif sfs_type == 'snp':
#
# 	maxiters = [5,20,100,None]
#
# 	initials = {
# 		"no_mig":[2,2,0.5],
# 		"sym_mig":[2,2,0.5,0.5],
# 		"sym_mig_size":[2,2,2,2,0.5,0.5,0.5],
# 		"anc_sym_mig":[2,2,0.5,0.5,0.5],
# 		"sec_contact_sym_mig":[2,2,0.5,0.5,0.5]
# 		}
#
# 	uppers = {
# 		"no_mig":[15,15,3],
# 		"sym_mig":[15,15,3,3],
# 		"sym_mig_size":[15,15,15,15,3,3,3],
# 		"anc_sym_mig":[15,15,3,3,3],
# 		"sec_contact_sym_mig":[15,15,3,3,3]
# 		}


param_labels = {
	"no_mig":"nu1, nu2, T",
	"sym_mig":"nu1, nu2, m, T",
	"CMG":"nu1, nu2, b1, b2, m, T",
	"anc_sym_mig":"nu1, nu2, m, T1, T2",
	"sec_contact_sym_mig":"nu1, nu2, m, T1, T2",
	"SCG":"nu1, nu2, b1, b2, m, Ts, Tsc",
	}

lowers = {
	"no_mig":[0.01,0.01,0.01],
	"sym_mig":[0.01,0.01,0.01,0.01],
	"CMG":[0.01,0.01,0.01,0.01,0.01,0.01],
	"anc_sym_mig":[0.01,0.01,0.01,0.01,0.01],
	"sec_contact_sym_mig":[0.01,0.01,0.01,0.01,0.01],
	"SCG":[0.01,0.01,0.01,0.01,0.01,0.01,0.01],
	}

#Set priors that are constant regardless of SFS type and model
pts_start = sum(ns)
pts = [pts_start, pts_start + 10, pts_start + 20]
if model == 'CMG' or model == 'SCG':
	reps = [20,20,30,20,5]
	folds = [3,2,1,1,0.5]
	maxiters = [1,5,20,50,None]
else:
	reps = [10,10,5,5]
	folds = [3,2,1,1]
	maxiters = [5,20,100,200]
func = getattr(Models_2D, model)
p_labels = param_labels[model]
in_params = initials[model]
upper = uppers[model]
lower = lowers[model]

#Create results_summary.txt file if this is the first run.
if not os.path.isfile("%s/results_summary/%s_%s_%s_%s_results_summary.txt" % (out_dir,taxa,is_folded,model,sfs_type)):
	summary_file = open("%s/results_summary/%s_%s_%s_%s_results_summary.txt" % (out_dir,taxa,is_folded,model,sfs_type), "a+")
	param_list = '\t'.join(map(str,param_labels[model].replace(" ","").split(",")))
	summary_file.write("sfs_type\trun_num\ttaxa\tmodel\tlog-likelihood\taic\ttheta\t%s\n" % (param_list))
	summary_file.close()

#Show grid points in case they need to be checked in the standard output file
print "\n\nGrid points: {}".format(pts)

print "\n\nBeginning optimization run {}".format(run_num)

#Set filename prefix
prefix = "%s/%s_%s_%s_run%s" % (out_dir,taxa,is_folded,sfs_type,run_num)

#Run optimization run with dadi_pipeline
if is_folded == "fold":
	params = Optimize_Functions.Optimize_Routine(fs, pts, prefix, model, func, len(reps), len(lower), fs_folded = True, param_labels = p_labels, reps = reps, folds = folds, in_params = in_params, in_upper = upper, in_lower = lower, maxiters = maxiters)
elif is_folded == "unfold":
	params = Optimize_Functions.Optimize_Routine(fs, pts, prefix, model, func, len(reps), len(lower), fs_folded = False, param_labels = p_labels, reps = reps, folds = folds, in_params = in_params, in_upper = upper, in_lower = lower, maxiters = maxiters)
else:
	print "\n\nSyntax incorrect for folding argument, please use 'fold' or 'unfold'"

#Print notable results and where they are saved to.
print "\n\nHighest likelihood of run: {}".format(params[0])
print "\nOptimized parameter values: {}".format(params[1])
print "\nAdding run results to analysis summary file: %s/results_summary/%s_%s_%s_%s_results_summary.txt" % (out_dir,taxa,is_folded,model,sfs_type)

#Once again make sure summary file is there.
if not os.path.isfile("%s/results_summary/%s_%s_%s_%s_results_summary.txt" % (out_dir,taxa,is_folded,model,sfs_type)):
	summary_file = open("%s/results_summary/%s_%s_%s_%s_results_summary.txt" % (out_dir,taxa,is_folded,model,sfs_type), "a+")
	param_list = '\t'.join(map(str,param_labels[model].replace(" ","").split(",")))
	summary_file.write("sfs_type\trun_num\ttaxa\tmodel\tlog-likelihood\taic\ttheta\t%s\n" % (param_list))
	summary_file.close()

#Write results to summary file.
summary_file = open("%s/results_summary/%s_%s_%s_%s_results_summary.txt" % (out_dir,taxa,is_folded,model,sfs_type),'a')
tab_params = '\t'.join(map(str,params[1]))
summary_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (sfs_type,run_num,taxa,model,params[0],params[2],params[3],tab_params))
summary_file.close()

print "\nRun complete..."

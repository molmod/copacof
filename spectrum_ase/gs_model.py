#!/bin/env python
# Siebe Vanlommel 2021 --..
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy.optimize import leastsq
import sys

shield = np.genfromtxt('./shieldings')

print('!annotated relative population only correct when $ntotalcarbons is set correctly!')
ntotalcarbons = 300

# some settings to change the way the distribution of the chemical shielding data is plotted
scatter_avgs = 1
boxplot_on = 1
showcaps = 0
showbox = 1
showfliers = 0
showmeans = 1
meanpointprops = dict(marker='D', markeredgecolor='black',
                      markerfacecolor='black', markersize=2)
show_difference_spectrum = 0 # plot difference between theoretical and experimental spectrum
annotate_labels = 0
boxplot_width = 0.01
boxplot_yLocations_scaler = 0.015

# plot parameters
# A0 = amplitude of spectrum, a0 and b0 the relation between shift and shielding and w0 the gamma parameter controlling overall width
# see main text for explanation
A0,a0,b0,w0 = 1.0, 1.0, 170.0, 1.0 # initial guess for fitting parameters in the gs model
Pars0 = [A0,a0,b0,w0]
fig,ax = plt.subplots(1, figsize = (10,8))

def gauss(freq, x, sigma):
   return 1.0/(sigma*np.sqrt(2.0*np.pi))*np.exp((-1.0/2.0)*((x-freq)/sigma)**2)

x_range = np.linspace(0,200,20000)
total_spectrum = np.array([0.0]*20000)
sigma = 1.0

nodes = np.genfromtxt('nodes_classified.dat', dtype='int')
linkers = np.genfromtxt('linkers_classified.dat', dtype='int')

node_colorfile = np.genfromtxt('nodes_colordef.dat') # change this file as you want the colors to be in different usecase
linker_colorfile = np.genfromtxt('linkers_colordef.dat')
all_colors = {}
for label, lr, lg, lb in node_colorfile:
	all_colors['n-{}'.format(int(label))] = (lr/256, lg/256, lb/256)
for label, lr, lg, lb in linker_colorfile:
	all_colors['l-{}'.format(int(label))] = (lr/256, lg/256, lb/256)

# make a spectrum by taking a gaussian per type
# with width standard deviation of shift data per class
# have to correct for weight of each class tho as not all classes are not necessarily equally populated
annotated = []
lineshapes_nodes = {}
lineshapes_linkers = {}
for ntype,ctype in enumerate(nodes):
	label = int(ctype[0])
	# determine the avg and the width of the lineshape
	type_shift = np.average(shield[ctype[1:]])
	type_width = np.std(shield[ctype[1:]])
	# add this as much as there are atoms in the type
	lineshapes_nodes[ntype] = [type_shift, type_width, len(ctype[1:])]

annotated = []
for ntype,ctype in enumerate(linkers):
	label = int(ctype[0])
	type_shift = np.average(shield[ctype[1:]])
	type_width = np.std(shield[ctype[1:]])
	lineshapes_linkers[ntype] = [type_shift, type_width, len(ctype[1:])]

####
# finally we can apply an overall scaling factor to the width of the lineshapes
# we'll fit this width as a parameter, together with the a and b coefficients
# in the shift-shielding relation, and the amplitude A
####
lineshapes = {}
for key,val in lineshapes_nodes.items():
	label = 'node-{}'.format(key)
	lineshapes[label] = val
for key,val in lineshapes_linkers.items():
	label = 'link-{}'.format(key)
	lineshapes[label] = val

# ls is the dictionary with all the lineshapes
# key :: node/link - *
# val :: [shift, intrinsic width, number in group]
def model(pars):
	comp = np.array([0.0]*len(x_ax))
	for label,shape in lineshapes.items():
		for count in range(shape[2]):
			comp += gauss(x_ax, shape[0] * pars[1] + pars[2], shape[1]*pars[3])

	return comp*pars[0]

# target is exp spectrum
exp = np.transpose(np.genfromtxt('spectrum.dat'))
x_ax = exp[0]
spec = exp[1]
spec -= np.min(spec)

# multiple options for choice of residual function

def residuals_mse(pars):
	ret = spec - np.array(model(pars))
	return ret

def residuals_rwp(pars):
    ret = spec - np.array(model(pars))
    weights = [1/speci if speci>0 else 0. for speci in spec]
    rwp = np.sqrt((weights*(ret)*(ret)).sum() /(weights*(spec*spec)).sum())
    return rwp

def residuals_rp(pars):
    ret = spec - np.array(model(pars))
    rp = np.sqrt(((ret)*(ret)).sum() /((spec*spec)).sum())
    return rp

def residuals_si(pars):
    p1 = spec
    p2 = np.array(model(pars))
    return (p1*p2).sum()/(np.sqrt((p1*p1).sum()) * np.sqrt((p2*p2).sum()))

# pick residual function
residual_function_choice = sys.argv[1]
if residual_function_choice == 'si':
	residual_function = residuals_si
elif residual_function_choice == 'mse':
	residual_function = residuals_mse
elif residual_function_choice == 'rwp':
	residual_function = residuals_rwp
elif residual_function_choice == 'rp':
	residual_function = residuals_rp

bounds = ([0.5e-2,0.8,150,0.1],[5e5,1.2,172,2.0]) # if desired, bounds may be given here to define the search space for the parameters
OptimizedModel = scipy.optimize.least_squares(fun = residual_function, x0 = np.array([A0,a0,b0,w0]), bounds = bounds)

Pars = OptimizedModel['x']
print('----')
print('optimal parameters: ')
print('overall amplitude: {:.3f}'.format(Pars[0]))
print('overall shift/shield: {:.3f}*shield + {:.3f}'.format(Pars[1],Pars[2]))
print('overall width scaling: {:.3f}'.format(Pars[3]))

best_model = model(Pars)

bplots = []
A_opt, a_opt, b_opt, w_opt = Pars[0], Pars[1], Pars[2], Pars[3]
for ntype,ctype in enumerate(nodes):
        label = int(ctype[0])
        type_shift = np.average(shield[ctype[1:]])
        type_width = np.std(shield[ctype[1:]])
        if scatter_avgs:
                ax.scatter(type_shift*a_opt + b_opt, ntype*boxplot_yLocations_scaler,
                    marker='D',alpha=0.6,color = all_colors['n-{}'.format(label)], s = 20)
        print('n-{} shift = {:.1f} ppm, width = {:.1f} ppm'.format(label,type_shift*a_opt + b_opt, type_width*w_opt))
        if annotate_labels:
                plt.annotate('n-{}({:.1f}%)'.format(label,len(ctype[1:]) / ntotalcarbons * 100), (type_shift*a_opt + b_opt, ntype*boxplot_yLocations_scaler),
                           #rotation=90, color=all_colors['n-{}'.format(label)])
                           color = 'black')
        if boxplot_on:
                ax.boxplot(shield[ctype[1:]]*a_opt + b_opt, vert=False, positions = [ntype*boxplot_yLocations_scaler], widths = [boxplot_width],
                    labels=['n-{}'.format(int(ctype[0]))], showcaps = showcaps,
                    showfliers = showfliers, showbox = showbox,
                    showmeans = showmeans, patch_artist=True,
                    boxprops=dict(facecolor=all_colors['n-{}'.format(label)]), meanprops = meanpointprops)

for ntype,ctype in enumerate(linkers):
        #current_col = next(bplot_cols)
        label = int(ctype[0])
        type_shift = np.average(shield[ctype[1:]])
        type_width = np.std(shield[ctype[1:]])
        if scatter_avgs:
                ax.scatter(type_shift*a_opt + b_opt, ntype*boxplot_yLocations_scaler,
                   marker='D',alpha=0.6,color = all_colors['l-{}'.format(label)], s = 20)
        print('l-{} shift = {:.1f} ppm, width = {:.1f} ppm'.format(label,type_shift*a_opt + b_opt, type_width*w_opt))
        if annotate_labels:
                plt.annotate('l-{}({:.1f}%)'.format(label,len(ctype[1:]) / ntotalcarbons * 100),
                  (type_shift*a_opt + b_opt, ntype*boxplot_yLocations_scaler),
                  #rotation=90, color=all_colors['l-{}'.format(label)])
                  color = 'black')
        if boxplot_on:
                ax.boxplot(shield[ctype[1:]]*a_opt + b_opt, vert=False, positions = [ntype*boxplot_yLocations_scaler],
                      widths = [boxplot_width], labels=['l-{}'.format(int(ctype[0]))], showcaps = showcaps,
                      showfliers = showfliers, showbox = showbox, showmeans = showmeans, patch_artist=True,
                      boxprops=dict(facecolor=all_colors['l-{}'.format(label)]), meanprops = meanpointprops)

#from sklearn.metrics import mean_squared_error as mse
print(np.trapz(spec, dx = 0.01))
#MSE = mse(spec, best_model) / np.trapz(spec, dx = 0.01)
MSE = np.mean((spec - best_model)**2) / np.trapz(spec, dx = 0.01)
print('MSE from exp: {} arbitrary units'.format(np.format_float_scientific(MSE, precision=2)))
ax.plot(x_ax, spec, label='exp')
ax.plot(x_ax, best_model, label='model, $g_s$', color = 'tab:red')
if show_difference_spectrum:
	ax.plot(x_ax, best_model - spec, label='diff', color='tab:green')
plt.legend()
ax.set_xlim((200.0,0.0))
ax.set_ylim((-0.05, np.max(spec) + 0.01))
ax.tick_params(labelsize=16)
ax.set_ylim((0.0 - np.max(spec)*0.025, np.max(spec)*1.025))
ax.set_yticks([])

plt.title('13C, residual_function = {} \n A={:.2f}, a={:.2f}, b={:.2f}, $\gamma$={:.2f} \n nMSE = {:.5f}'.format(residual_function_choice, A_opt, a_opt, b_opt, w_opt, MSE))

plt.show()

import getdist.plots as gplot
import os


analysis_settings = {'ignore_rows': u'0.3'}
roots = ['pk15_VVE','pk15_LSS_VVE']
#roots= ['Planck_LSS_bigprior_varbin']
params = ['H0','void_qV0','omegam','void_z0']
colors = ['#8E001C','#FFB300','navy']
labels = ['Planck','Planck+LSS']

param_3d = None
g=gplot.getSubplotPlotter(chain_dir=r'./chains',analysis_settings=analysis_settings)

g.triangle_plot(roots, params, contour_colors=colors, legend_colors=colors, legend_labels=labels, plot_3d_with_param=param_3d, filled=True, shaded=False)

g.export('results_plots/trivar_VVE.pdf')

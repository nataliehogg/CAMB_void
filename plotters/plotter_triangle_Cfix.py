import getdist.plots as gplot
import os


analysis_settings = {'ignore_rows': u'0.3'}
roots = ['pk15_Cfix','pk15_LSS_Cfix']
#roots= ['pk15_LSS_Cfix']
params = ['H0','void_qV0','omegam']
colors = ['#8E001C','#FFB300','navy']
#colors = ['#FFB300']
labels = ['Planck','Planck+LSS','Planck bigstep']
#labels = ['Planck+LSS']

param_3d = None
g=gplot.getSubplotPlotter(chain_dir=r'./chains',analysis_settings=analysis_settings)

g.triangle_plot(roots, params, contour_colors=colors, legend_colors=colors, legend_labels=labels, plot_3d_with_param=param_3d, filled=[True,True,False], shaded=False)

g.export('results_plots/trivar_Cfix.pdf')

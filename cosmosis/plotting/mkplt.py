#Figure 4 has our fiducial chain, non-tomographic constraints, 
#pseudo Cls, Cl bandpowers, and Planck
from cosmosis.postprocessing import plots
from cosmosis.postprocessing import lazy_pylab as pylab
from cosmosis.postprocessing import statistics
from cosmosis.plotting.kde import KDE

from cosmosis.postprocessing.elements import PostProcessorElement
from cosmosis.postprocessing.elements import MCMCPostProcessorElement, MultinestPostProcessorElement, WeightedMCMCPostProcessorElement
from cosmosis.postprocessing.elements import Loadable
from cosmosis.postprocessing.outputs import PostprocessPlot
from cosmosis.postprocessing.utils import std_weight, mean_weight

import numpy as np
import scipy

import matplotlib.colors
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['font.size']=18
matplotlib.rcParams['legend.fontsize']=50
matplotlib.rcParams['xtick.major.size'] = 10.0
matplotlib.rcParams['ytick.major.size'] = 10.0

class plot2D(plots.MetropolisHastingsPlots2D):
    contours=[]
    proxies=[]
    def __init__(self, *args, **kwargs):
        super(plot2D, self).__init__(*args, **kwargs)
        self.colors=['purple', "slateblue", "forestgreen", "pink"]
        self.linestyles=[':','--','-','-']
        self.labels=[ "$\Delta m = 0.05$", "$\Delta m = 0.02$", "$\Delta m = 0.01$", "Line 4"]
        self.axis=[0.2,0.45, 0.65,1.0] #[0.233,0.45,0.65,1.1]
        self.fill_list=[True,True,True,True]*10
        self.line_list=[True,True,True,True]*10
        self.opaque_list=[True,True,True,True]*10
        self.linewidths=[2.5]*len(self.colors)
        self.alphas=[0.4,0.4,0.4,0.4]*len(self.colors)
        self.imshow=[False]*10
        self.linecolors=[None]*10
        self.opaque_centre=[False]*10
        self.fill_colors=[None,None]*10
        self.plot_points=[(0.3156,0.82)]
        #self.proxies=[]

    def run(self):
        name2="cosmological_parameters--sigma8_input"
        name1="cosmological_parameters--omega_m"
        if "fisher" not in self.source.name:
            filename = self.make_2d_plot(name1, name2)
        else:
            filename=None  
        return filename

    def smooth_likelihood(self, x, y):
        n = self.options.get("n_kde", 100)
        factor = self.options.get("factor_kde", 2.0)
        kde = KDE([x,y], factor=factor)
        x_range = (x.min(), x.max())
        y_range = (y.min(), y.max())
        (x_axis, y_axis), like = kde.grid_evaluate(n, [x_range, y_range])
        return n, x_axis, y_axis, like
    def make_2d_plot(self, name1, name2):
          #Get the data
        x = self.reduced_col(name1)
        y = self.reduced_col(name2)


        if x.max()-x.min()==0 or y.max()-y.min()==0:
            return
        print "  (making %s vs %s)" % (name1, name2)


        #Interpolate using KDE
        try:
            n, x_axis, y_axis, like = self.smooth_likelihood(x, y)
        except np.linalg.LinAlgError:
            print "  -- these two parameters have singular covariance - probably a linear relation"
            print "Not making a 2D plot of them"
            return []

        figure,filename = self.figure("2D", name1, name2)

          #Choose levels at which to plot contours
        
        contour1=1-0.68
        contour2=1-0.95
        try: 
            level1, level2, total_mass = self._find_contours(like, x, y, n, x_axis[0], x_axis[-1], y_axis[0], y_axis[-1], contour1, contour2)
        except:
            import pdb ; pdb.set_trace()
        level0 = 1.1
        levels = [level2, level1, level0]
        #Create the figure
        fig,filename = self.figure("2D", name1, name2)
        pylab.figure(fig.number)

        #Make the plot
        #pylab.figure(figure.number, figsize=(8,7))
        keywords = self.keywords_2d()
        fill = self.options.get("fill", True)
        if self.fill_list!=None:
            fill = self.fill_list[self.source.index]
        imshow = self.imshow[self.source.index]
        #plot_points = self.options.get("plot_points", False)
        color = self.colors[self.source.index]
        linestyle = self.linestyles[self.source.index]
        if self.linecolors[self.source.index] is not None:
            linecolor=self.linecolors[self.source.index]
        else:
            linecolor=color
        print self.source.index,linecolor
        
        if self.opaque_list[self.source.index]:
            pylab.contourf(x_axis, y_axis, like.T, [level2,level0], colors=['white'])
        print "fill",self.source.index, fill

        if imshow:
            pylab.imshow(like.T, extent=(x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]), aspect='auto', origin='lower', cmap=cmap_white_to_color(color))
        elif fill:
            if self.fill_colors[self.source.index] is not None:
                c=self.fill_colors[self.source.index]
                pylab.contourf(x_axis, y_axis, like.T, [level2,level0], colors=[c[0]])
                pylab.contourf(x_axis, y_axis, like.T, [level1,level0], colors=[c[1]])
            else:
                pylab.contourf(x_axis, y_axis, like.T, [level2,level0], colors=[color], alpha=self.alphas[self.source.index])
                if self.opaque_centre[self.source.index]:
                    pylab.contourf(x_axis, y_axis, like.T, [level1,level0], colors=[color])
                else:
                    pylab.contourf(x_axis, y_axis, like.T, [level1,level0], colors=[color], alpha=self.alphas[self.source.index])            
        if (not fill) or self.line_list[self.source.index]: 
            print self.source.index,linecolor
            pylab.contour(x_axis, y_axis, like.T, [level2,level1], colors=[linecolor], linestyles=linestyle, linewidths=self.linewidths[self.source.index])
            #if self.labels is not None:
            #    self.proxies.append(pylab.plot([],[], color=color, linestyle=linestyle, linewidth=self.linewidths[self.source.index])[0])
        if self.labels is not None:
            if self.fill_colors[self.source.index] is not None:
                print self.fill_colors[self.source.index]
            elif fill:
                print self.source.index,self.line_list[self.source.index]
                if self.line_list[self.source.index] is not None:                    
                    conv=matplotlib.colors.ColorConverter()
                    edgecolor=linecolor
                    facecolor=conv.to_rgba(color,alpha=1-(1-self.alphas[self.source.index])**2)
                    self.proxies.append(pylab.Rectangle((0,0),1,1,facecolor=facecolor,edgecolor=edgecolor))#,alpha=(1-(1-self.alphas[self.source.index])**2)))
                else:
                    self.proxies.append(pylab.Rectangle((0,0),1,1,fc=color,edgecolor=color,alpha=(1-(1-self.alphas[self.source.index])**2)))
            else:
                self.proxies.append(pylab.plot([],[], color=color, linestyle=linestyle, linewidth=self.linewidths[self.source.index])[0])

        if self.plot_points is not None:
            point_markers=['x','+','^']
            for p,m in zip(self.plot_points,point_markers):
                print p,m
                pylab.plot(p[0], p[1], m, color='k', markersize=10, markerfacecolor='none')

        #Do the labels
        print self.proxies,self.labels
        if self.labels is not None:
            leg=pylab.legend(self.proxies,self.labels,loc="upper right")
            leg.get_frame().set_alpha(0) # this will make the box totally transparent
            leg.get_frame().set_edgecolor('white') # this will make the edges of the 
        pylab.xlabel(self.latex(name1))
        pylab.ylabel(self.latex(name2))
        pylab.tight_layout()
        if self.axis is not None:
            pylab.axis(self.axis)
        return filename     

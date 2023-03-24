"""
Applies the SNHA on the Vibrio - time-lagged environmental dataset and plots
the network (Figure 6).
"""

from Snha import Snha
import pandas as pd

if __name__ == "__main__":
    
    # load the data
    fpath = "N:/data/merged_vibrio_env/"
    df = pd.read_csv(fpath + "vibrio_data_merged_sh_mv_timelag.csv", 
                     sep = "\t")
    
    s = Snha(df)  # initialize the St Nichlaus object with the data

    # compute the correlation with Spearman rank correlation
    s.compCorr(method = "spearman")
    
    # run the algorithm with bootstrap
    s.stNicholasAlg(alpha = 0.1, # correlation coefficient cut off;default: 0.1
                    bt = True, # Bootstrap
                    n = 100, #  number of bootstrap iterations; default: 100
                    lbd = 0.08, # fraction of all iterations to accept an edge
                                # as a prediction; default: 0.5 (e.g. 50/100 
                                # iterations an edge need to be found, to be 
                                # considered as a predicted edge at lbd=0.5 
                                # and n=100)
                    method = "spearman",
                    )
    
    # # optional: confirm that selected lbd = 0.08 in combination with n = 100
    # # returns only significant edges:
    # from scipy.stats import binom_test
    # binom_test(8,100,0.035)
    # # 0.035 is the likely the highest probability of success 
    # # (i.e. chance of detecting false positives. Depends on the number of 
    # # input varaibles)
    
    #%% plot the graph
    import igraph as ig
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib import cm
    import cmocean
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import random

    class _cmapHelper:
      def __init__(self, cmap = cmocean.cm.balance, vmin = -1, vmax = 1):
        self.cmap = cmap
        self.norm = mpl.colors.Normalize(vmin=-1, vmax=1)
        self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)
      def get_rgb(self, val):
        return self.scalarMap.to_rgba(val)

    def _insertCbar(ax, cmap, label = "Spearman R", orientation = "vertical"):
        
        if orientation == "horizontal":
            cax = inset_axes(ax,
                               width="50%",  
                               height="5%",
                               loc='lower center',
                               borderpad=-1
                               )
            
        elif orientation == "vertical":
            cax = inset_axes(ax,
                               width="5%",  
                               height="50%",
                               loc='center right',
                               borderpad=-5
                               )
            
        else:
            raise ValueError("Parameter 'orientation' must be one of ['horizontal'"
                             ", 'vertical']. Is {} instead.".format(orientation))
                             
        plt.colorbar(cmap.scalarMap, 
                     cax = cax, 
                     label = label, 
                     orientation = orientation)
    
    # define random seed (defines how the variables are aranged)
    random.seed(2)#7

    labels = ['Vibrio\nqty', #'quantity', 
               'Chl\nmean\n6-15', # different from max r_s (15\n15)
               'Chl\ntrend\n3-5', 
               'O2\nmean\n0-0',
               'O2\ntrend\n5-16', 
               'NO3\nmean\n7-25', 
               'NO3\ntrend\n0-24', 
               'NH4\nmean\n28-28',
               'NH4\ntrend\n9-10', 
               'PO4\nmean\n22-22',
               'PO4\ntrend\n4-9', # different from max r_s (20\n21)
               'SST\nmean\n0-11',
               'SST\ntrend\n0-30', # different from max r_s (24\n25)
               'SST\n180',
               'SSS\nmean\n0-6',
               'SSS\ntrend\n11-26',
               'Irrad\nmean\n27-28',
               'Irrad\ntrend\n18-28',
               'SAT\nmean\n0-16',
               'SAT\ntrend\n7-29',
               'Prec\nmean\n6-18',
               'Prec\ntrend\n5-8', 
               'WS\nmean\n2-17',
               'WS\ntrend\n4-8',
              ]    

    # create the matplotlib figure
    fig, ax = plt.subplots(figsize=(8, 8))

    # create the Graph from the adjacency matrix
    G = ig.Graph.Adjacency(s.graphPred, mode="max")
    
    # hide the axis
    plt.axis('off')

    # define the layout of the graph
    layout = G.layout_davidson_harel()

    # get the correlation matrix as a df
    corr_df = s.corr
    # initialize the cmapHelper class
    cmap = _cmapHelper()
    
    # define the visual style of the iGraph
    visual_style = {}
    visual_style["vertex_color"] = "Gray"
    visual_style["vertex_label"] = labels
    # visual_style["vertex_size"] = 0.8
    visual_style["vertex_size"] = 1.7
    visual_style["vertex_label_size"] = 7
    edge_corrs = [corr_df.iloc[edge] for edge in G.get_edgelist()]
    
    # get RGBs for each edge based on the cmap in the cmapHelper-object
    visual_style["edge_color"] = [cmap.get_rgb(x) for x in edge_corrs]
    visual_style["edge_width"] = [4*abs(x) for x in edge_corrs]
    
    # plot the graph
    ig.plot(
        G,
        layout=layout,
        bbox = (2000, 2000), # defines the size of the plot in pixels!
        target = ax,
        **visual_style
    )
    
    # insert the colorbar
    _insertCbar(ax, cmap, label = "$r_s$", orientation = "horizontal")
    
    plt.savefig("N:/plots/snha_bootstrap_100_it.pdf",
                dpi = 300,
                bbox_inches = "tight")
    
    plt.show()
    

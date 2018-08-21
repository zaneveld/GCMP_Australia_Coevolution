import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sys import argv


def make_heatmap(data,outfile,title=None, column_normalization=None,cell_value='R2_adj',cmap='vlag',annot=True):
    """Generate and save a heatmap

    data -- a Pandas datafram with columns called 'category','Compartment' and whatever is
    specified in cell_value (typically 'R2' or 'R2_adj')

    outfile -- where to save the figure image

    title -- title of plot (string)

    cell_value -- column to use as data when making the clustermap

    column_normalization -- if 1, use z-score normalization on columns (asking which row best
    explains column). if 0 normalize by rows. Use None to avoid normalization.

    cmap - colormap to use, e.g. 'viridis','plasma','BuGn','mako', etc.

    annot - add values to cells (values will be normalized if column_normalization is true
    """    
    curr_data = data.loc[data['bdiv_metric']==metric,['category','Compartment',cell_value]]
    pivoted_data = curr_data.pivot(index='Compartment',columns='category',values=cell_value) 
    #Add row colors
    pivoted_data = pivoted_data.reindex(index=['mucus','tissue','skeleton'])    
    
    #This version is normalized by column z-score
    #addressing "what compartment most clusters each fcator"    
    g = sns.clustermap(pivoted_data,annot=annot,metric="correlation",cmap=cmap,\
     figsize=(24,8),row_cluster=True,robust=False,z_score=column_normalization)

 
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90,fontsize=8)
    plt.subplots_adjust(bottom=0.38)
    plt.suptitle(title)
    #plt.show()
    g.savefig(outfile,format='pdf')
    plt.clf()

if __name__ == "__main__":
    infile = argv[1]
    data = pd.read_table(infile)
    cell_value = 'R2_adj'
    #cmap = 'RdPu'
    cmap = 'Wistia'
    sns.set_context('paper')
    for column_normalization in [None,1]:
        
        col_norm_method = 'raw'
        if column_normalization == 1:
            col_norm_method = 'z_score_normalized'
        
        
        for metric in ['unweighted_unifrac','weighted_unifrac','bray_curtis']:
            for annot in [True,False]:
                title = "Adonis %s %s %s" %(cell_value,metric,str(col_norm_method))
                outfile = '%s_heatmap_%s_%s_%s_annot_%s.pdf'%(infile,cell_value,metric,col_norm_method,str(annot))
                make_heatmap(data,outfile,title,column_normalization,cell_value,annot=annot,cmap=cmap)
                 

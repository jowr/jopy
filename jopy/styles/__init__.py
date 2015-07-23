
import matplotlib.pyplot as plt 

from .plots import Figure

def get_figure(orientation='landscape',width=110,fig=None,axs=False):
    """Creates a figure with some initial properties
    
    The object can be customised with the parameters. But since it is an 
    object, it can also be modified later.
    
    Parameters
    ----------
    orientation : str
        either landscape or portrait
    width : float
        width in mm, used for scaling the default A4 format
        if width is less than 10, use it as a factor to apply to A4
    fig : matplotlib.figure.Figure, jopy.plots.Figure
        The figure object to handle, None creates a new figure
    axs : boolean
        True or False - Should an axis object be added to the figure?
    
    Returns
    -------
    jopy.plots.Figure, matplotlib.figure.Figure
        The figure object
    
    """

    if fig is None: fig = plt.figure()
    
    fig.__class__ = Figure 

    sideA = 297. # height of A4
    sideB = 210. # width of A4
    mm_to_inch = 3.93700787401575/100.0 # factor mm to inch

    if width < 0: raise ValueError("size cannot be less than zero.")

    width *= mm_to_inch
    sideA *= mm_to_inch
    sideB *= mm_to_inch

    if orientation=='landscape':
        if width < 10*mm_to_inch: width *= sideA
        scale  = width/sideA
        width  = sideA*scale #=width
        height = sideB*scale
    elif orientation=='portrait':
        if width < 10*mm_to_inch: width *= sideB
        scale  = width/sideB
        width  = sideB*scale #=width
        height = sideA*scale
    else:
        raise ValueError("Unknown orientation")
    fig.set_size_inches(width,height)
    if axs: fig.add_subplot(111)
    return fig


def plot_axis(data,kind,ax=None):
    if ax is None: ax = plt.gca()
    




if __name__ == "__main__":
    from jopy.utils import module_class_dict
    import jopy.styles.mplib as mpl
    dic = module_class_dict(mpl)
    for i in dic:
        line_fig,map_fig = dic[i]()._show_info()
        line_fig.savefig(i+"_lines.pdf")
        map_fig.savefig(i+"_maps.pdf")

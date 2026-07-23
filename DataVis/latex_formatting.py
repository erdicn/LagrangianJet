# Use the below only once the plots are clearly done so it is ready for presentation style
# This cell was taken from gemini and modified accordingly
# import matplotlib as mpl

# mpl.use("pgf")
# latex_textwidth_inches = 5.9 
# plt.rcParams.update({
#     "pgf.texsystem": "pdflatex",   # or 'xelatex' / 'lualatex'
#     "text.usetex": True,           # use LaTeX to write all text
#     "font.family": "serif",        # use serif (matching LaTeX default)
#     "pgf.rcfonts": False,          # don't setup fonts from rc parameters
#     "font.size": 11,               # MATCH YOUR LATEX DOCUMENT FONT SIZE
#     "axes.labelsize": 11,          # Match LaTeX font size
#     "legend.fontsize": 9,          # Slightly smaller for legends
#     "xtick.labelsize": 9,          # Slightly smaller for ticks
#     "ytick.labelsize": 9,
#     # CRITICAL: Set figure width to match LaTeX text width!
#     "figure.figsize": (latex_textwidth_inches, 4) # (width, height)
# })

import matplotlib.pyplot as plt
import re

def configureLatexPlot(context="twocolumn", aspect_ratio=None, caption_spacer=0.7, column_size = None, bigger_font_by = 0):
    """
    Configures Matplotlib rcParams for perfect LaTeX exports.
    
    Parameters:
    context (str): 'article' (full width), 'twocolumn' (half width), or 'beamer' (slides)
    """
    
    # 1. Define dimensions and fonts based on standard LaTeX templates
    if context == "twocolumn":
        # Typical IEEE/ACM two-column format
        # \columnwidth is usually ~252pt (3.5 inches)
        fig_width = 3.5  
        fig_height = 2.5 
        base_font_size = 11 # Two-column papers usually use 9pt or 10pt base fonts
        
    elif context == "beamer":
        # Standard Beamer slide (4:3 aspect ratio)
        # \textwidth is usually ~307pt (4.25 inches)
        latex_width_pt  = 398.3386
        latex_height_pt = 233.74188
        latex_width_inch  = latex_width_pt/72.27
        latex_height_inch = latex_height_pt/72.27
        if column_size is not None:
            latex_width_inch *= column_size
        fig_width  = latex_width_inch#4.25 
        fig_height = latex_height_inch*caption_spacer # multiplier so i have space for the caption#2.75 
        base_font_size = 10 # Beamer fonts look big on screen, but base is 11pt
    
    elif context == "halfbeamer":
        # Standard Beamer slide (4:3 aspect ratio)
        # \textwidth is usually ~307pt (4.25 inches)
        latex_width_pt  = 398.3386/2
        latex_height_pt = 233.74188
        latex_width_inch  = latex_width_pt/72.27
        if column_size is not None:
            latex_width_inch *= column_size
        latex_height_inch = latex_height_pt/72.27
        fig_width  = latex_width_inch#4.25 
        fig_height = latex_height_inch*caption_spacer# multiplier so i have space for the caption #2.75 
        base_font_size = 11 # Beamer fonts look big on screen, but base is 11pt    
        # plt.rcParams.update({
        #     "legend.handlelength":1,
        #     "handletextpad":0.4,
        #     "borderpad":0.4
        #     # plt.legend(
        #     #     loc="upper left",
        #     #     fontsize=8, 
        #     #     handlelength=1.0, 
        #     #     handletextpad=0.4, 
        #     #     borderpad=0.4
        #     # )
        # })
    elif context == "article":
        # Standard single column (e.g., typical report/thesis)
        # \textwidth is usually ~426pt (5.9 inches)
        fig_width = 5.9  
        fig_height = 4.0 
        base_font_size = 11
        
    else:
        raise ValueError("Context must be 'article', 'twocolumn', 'halfbeamer', or 'beamer'")

    if aspect_ratio is not None:
        fig_height = fig_width/aspect_ratio
    
    # 2. Apply the settings globally
    base_font_size += bigger_font_by
    
    plt.rcParams.update({
        # "text.usetex": True,
        # "font.family": "serif",
        # "font.serif": ["Computer Modern Roman"], 
        
        # Scaling everything based on the base_font_size
        "font.size": base_font_size,
        "axes.titlesize": base_font_size,
        "axes.labelsize": base_font_size,
        "legend.fontsize": base_font_size - 2,   # Slightly smaller for space
        "xtick.labelsize": base_font_size - 2,   # Smaller tick labels
        "ytick.labelsize": base_font_size - 2,
        
        # Set exact figure dimensions
        "figure.figsize": (fig_width, fig_height),
        
        # Automatically adjust padding so labels don't get cut off
        # "figure.autolayout": True 
        "figure.constrained_layout.use": True
    })
    
    if re.search(r"beamer", context, re.IGNORECASE):
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "sans-serif",
            "font.sans-serif": ["Computer Modern Sans Serif"],
            # This makes math equations sans-serif, matching Beamer's default math style
            "text.latex.preamble": r"\usepackage{sansmath} \sansmath"
        })
    else :
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"], 
            "text.latex.preamble": r"\usepackage{amsmath}"
            
        }) 
    
    print(f"Plot configured for {context}: {fig_width}\"x{fig_height}\" at {base_font_size}pt.")

def backToDefaultMatplotlib():
    plt.rcdefaults()
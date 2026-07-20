import numpy as np
rho_water = 998
def volumeSphere(r):
    return 4/3*np.pi*r**3
def massSphere(r, rho = rho_water):
    return rho*volumeSphere(r)
def massToRadius(m, rho=rho_water):
    return np.cbrt(m/rho * 3/4 /np.pi)

def getProbDenisty(raw_data, bin_nb = 300):
    # counts, bin_edges = np.histogram(raw_data, bins=bin_nb, density=True)
    probabilities, bin_edges = np.histogram(raw_data, bins=bin_nb, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2 
    # probabilities = counts / counts.sum() # This array will be bound between 0 and 1
    return probabilities, bin_centers, bin_edges

def getProbDenistyByVolume(diameter_data, bin_nb = 300):
    # counts, bin_edges = np.histogram(diameter_data,weights=volumeSphere(diameter_data/2), bins=bin_nb, density = True)
    probabilities, bin_edges = np.histogram(diameter_data,weights=volumeSphere(diameter_data/2), bins=bin_nb, density = True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2 
    # probabilities = counts / counts.sum() # This array will be bound between 0 and 1
    return probabilities, bin_centers, bin_edges

def probDistByVolToProbDist(prob_dist_by_vol, d):
    f = prob_dist_by_vol/volumeSphere(d/2)
    return f/np.trapezoid(f, d)
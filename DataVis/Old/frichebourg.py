import numpy as np

def frichebourgF(Mp, omega1 = 2.33811):
    return 2*np.sqrt(np.pi)*pow(Mp,(5/2))*np.exp(-Mp*omega1 - (pow(Mp,3))/12)

def frichebourgFSmallM(Mp):
    return 1/(np.sqrt(np.pi*Mp))

def frichebourgAnaRho(m, t, Mp, F):
    return m*m/(pow(t,(4/3)))*F(Mp)

def frichebourgAnaRhoFit(M, m, t):
    Mp = M/(t**(2/3))
    return m*m/(pow(t,(4/3)))*frichebourgF(Mp)


"""#EXAMPLE USAGE
guess = [0.7, 10.0]
# 2. Add physical boundaries (m and t must be positive, and maybe t shouldn't exceed 50)
lower_bounds = [0.001, 0.001]
upper_bounds = [10.0, 50.0]

popt, pcov = curve_fit(
    frichebourgAnaRhoFit, 
    bin_centers[nn], 
    probabilities[nn], 
    p0=guess, 
    bounds=(lower_bounds, upper_bounds)
)

print(f"Optimized from your guess: m={popt[0]:.3f}, t={popt[1]:.3f}")
m_opt, t_opt = popt
plt.plot(bin_centers, frichebourgAnaRhoFit(bin_centers, m_opt, t_opt), label = "frichebourg fit all data")
"""


# For gamma fit 
"""

shape, loc, scale = gamma.fit(r_data)
gamma_pdf = gamma.pdf(bin_centers, shape, loc=loc, scale=scale)
scaled_gamma = gamma_pdf * bin_width
plt.plot(bin_centers, scaled_gamma, label=f'Gamma Fit (shape={shape:.2f})', linestyle='--')

"""
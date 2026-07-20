from data_analysis import *
from data_extraction import getDataFromFolder
class MyProbs:
    def printSurface(self):
        to_plot = self
        # print("area under actual diameter prob =", np.trapezoid( to_plot.d_actual_probability,     to_plot.d_prob_dist_bin_centers ))
        # print("area under m actual prob        =", np.trapezoid( to_plot.m_actual_probability,     to_plot.m_prob_dist_bin_centers))
        print("area under diameter prob dist   =", np.trapezoid( to_plot.d_prob_density         ,     to_plot.d_prob_dist_bin_centers ))
        print("area under m distribution       =", np.trapezoid( to_plot.m_prob_density,              to_plot.m_prob_dist_bin_centers))
        if self.d_vol_weighted_ptob_dist is not None:
            print("area under d prob*volume        =", np.trapezoid( to_plot.d_vol_weighted_ptob_dist, to_plot.d_prob_weighted_dist_bin_centers))
        
    def fromConcentration(self, diameters, concentration, r = 1.047, center = "geo"):
        #geometric center 
        if center == "geo":
            di = diameters/np.sqrt(r)
            dip1 = diameters*np.sqrt(r) 
        elif center == "ari": # arithmetic
            di = diameters*(2/(1+r))
            dip1 = diameters*2*r/(1+r) 
        else:
            raise ValueError("center type must be 'geo' or 'ari'")
            print("diameters are taken as bin centers")
            di   = diameters[:-1]
            dip1 = diameters[1:]
            # probably throws an error because we need to rescale diameters to be compatible with the same size
        delta_i = dip1-di
        
        n_d_i_delta_i = concentration
        self.concentration = concentration
        n_d_i = n_d_i_delta_i/delta_i 
        norm_factor = np.sum(n_d_i*delta_i)
        
        self.d_prob_density = n_d_i/norm_factor
        self.d_prob_dist_bin_centers = diameters
        # g    mg  µg
        # 1_000_000_000 / (10^6)^3 = 10^9/10^18 = 10^-9
        m_lower  = massSphere(di/2       , rho = rho_water*1e-9)
        m_upper  = massSphere(dip1/2     , rho = rho_water*1e-9)
        m_center = massSphere(diameters/2, rho = rho_water*1e-9)
        
        delta_m_i = m_upper - m_lower
        
        mass_in_bin_i = concentration * m_center # not exact but good enough ? 
        
        mass_density_i = mass_in_bin_i / delta_m_i
        
        m_norm_factor = np.sum(mass_density_i * delta_m_i)
        pdf_m = mass_density_i / m_norm_factor
        relative_mass_i = concentration * (diameters**3)
        mass_density_i = relative_mass_i / delta_m_i         
        m_norm_factor = np.sum(mass_density_i * delta_m_i) 
        
        self.m_prob_density = mass_density_i / m_norm_factor
        self.m_prob_dist_bin_centers = m_center

        d_area = np.sum(self.d_prob_density * delta_i)
        print("Diameter Area:", d_area)
        m_area = np.sum(self.m_prob_density * delta_m_i)
        print("Mass Area:", m_area)

        return 
    def __init__(self, diameters, label, n_bins = 500, experimental_prob_vol = None, 
                 concentration = None, velocities=None, linestyle = "-", experimental_vel_distrib_sat_main_x_y=None): # put experimental_vel_distrib_sat_main_x_y = ((satelitex, satelitey), (mainx, mainy))
        self.label = label
        self.d = diameters
        self.m = massSphere(self.d/2, rho = rho_water)
        self.linestyle = linestyle
        print(label)
        
        if velocities is not None:
            self.v = velocities
            print("NEED THE DIAMETERS IN METERS")
            self.v_main_prob_density, self.v_main_prob_dist_bin_centers, _         = getProbDenisty(self.v[self.d>75e-6], n_bins)
            self.v_satelite_prob_density, self.v_satelite_prob_dist_bin_centers, _ = getProbDenisty(self.v[self.d<75e-6], int(n_bins/2))
        
        if concentration is not None:
            self.fromConcentration(diameters, concentration)
            self.d_vol_weighted_ptob_dist        = None
            self.d_prob_weighted_dist_bin_centers = None
            return
        
        if experimental_prob_vol is not None:
            self.d_vol_weighted_ptob_dist = experimental_prob_vol 
            self.d_prob_dist_bin_centers  = diameters
            self.d_prob_weighted_dist_bin_centers = diameters
            self.d_prob_density =  probDistByVolToProbDist(experimental_prob_vol, diameters) 

            dd_bin_widths = np.gradient(self.d)
            self.dd_bin_widths =dd_bin_widths
            self.d_actual_probability = self.d_prob_density * dd_bin_widths

            self.m_prob_dist_bin_centers = self.m 
            dm_dd = rho_water * (np.pi / 2.0) * (self.d)**2

            # dx en fonction de d 
            self.m_prob_density = self.d_prob_density / dm_dd 
            dm_bin_widths = np.gradient(self.m)
            self.m_actual_probability = self.m_prob_density * dm_bin_widths
            
            if experimental_vel_distrib_sat_main_x_y is not None:
                # sat_x, sat_y   = experimental_vel_distrib_sat_main_x_y[0]
                # main_x, main_y = experimental_vel_distrib_sat_main_x_y[1]
                self.v_satelite_prob_dist_bin_centers, self.v_satelite_prob_density = experimental_vel_distrib_sat_main_x_y[0]
                self.v_main_prob_dist_bin_centers    , self.v_main_prob_density     = experimental_vel_distrib_sat_main_x_y[1] 
            
            return

        self.n_bins = n_bins
        
        self.m_prob_density, self.m_prob_dist_bin_centers, _ = getProbDenisty(self.m, n_bins)
        # self.m_vol_weighted_ptob_dist, self.m_prob_dist_bin_centers, _ = getProbDenistyByVolume(self.m, n_bins)
        
        self.d_prob_density, self.d_prob_dist_bin_centers, _ = getProbDenisty(self.d, n_bins)
        self.d_vol_weighted_ptob_dist, self.d_prob_weighted_dist_bin_centers, _ = getProbDenistyByVolume(self.d, n_bins)
        
        self.printSurface()
        
    @classmethod
    def initFromFolder(cls, folder_path, label, n_bins=500, linestyle='-', diameters_to_si=1e-6, vel_to_si=1):
        return cls(diameters  = diameters_to_si*getDataFromFolder(get_from_dir=folder_path, file_name_regex="radius*.dat")*2,
                   velocities = vel_to_si      *getDataFromFolder(get_from_dir=folder_path, file_name_regex="velocity*.dat"),
                   label=label)
        
    @classmethod
    def initFromFolderList(cls, l_folder_paths, label, n_bins=500, linestyle='-', diameters_to_si=1e-6, vel_to_si=1):
        diameters = []
        velocities= [] 
        for f in l_folder_paths:
            diameters .append(diameters_to_si*2*getDataFromFolder(get_from_dir=f, file_name_regex="radius*.dat"))
            velocities.append(vel_to_si      *  getDataFromFolder(get_from_dir=f, file_name_regex="velocity*.dat"))
        diameters = np.vstack(diameters)
        velocities = np.vstack(velocities)
        return cls(diameters  = diameters, velocities = velocities, label=label)
    
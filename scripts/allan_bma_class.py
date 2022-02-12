class bmaValidation:
    """This class provides a set of metrics and plots to validate the bma code
    """
    def __init__(self):
        print('Welcome to BMA Validation')
    
    def add_model(self, name, mass, mass_true, abs_mag, abs_mag_true, redshift):
        self.model_name = name
        self.mass = np.array(mass)
        self.mass_true = np.array(mass_true)
        self.Mr = np.array(abs_mag)
        self.Mr_true = np.array(abs_mag_true)
        self.z = np.array(redshift)
        
        # new variables
        self.res_mass = self.mass-self.mass_true
        self.res_Mr = self.Mr-self.Mr_true
                
    def plot_residual_mass(self,name,mask=None):
        if mask is None: mask = np.argsort(self.mass)
        ylabel = res_mass_label
        # put the plot function you made here
        pass
    
    def plot_resiudal_absolute_mag(self,mask=None):
        if mask is None: mask = np.argsort(self.Mr)
        ylabel = res_abs_mag_label
        # put the plot function you made here
        pass
    
    def plot_identity_mass_redshift(self,mask=None):
        if mask is None: mask = np.argsort(self.mass)
        tree_variable_plot(self.mass_true[mask],self.mass[mask],self.z[mask])
        plt.ylabel(res_mass_label)
        plt.xlabel(mass_true_label)
        pass

    def plot_identity_mass_chisqr(self,mask=None):
        if mask is None: mask = np.argsort(self.mass)
        tree_variable_plot(self.mass_true[mask],self.mass[mask],self.chisqr[mask],zlabel='Log(chisqr)')
        plt.ylabel(res_mass_label)
        plt.xlabel(mass_true_label)
        pass

# auxialiary variables and functions
mass_label = r'Log($M_{\star}^{BMA}$)'
mass_true_label = r'Log($M_{\star}^{COSMOS}$)'
res_mass_label = r'Log $\left(M_{\star}^{COSMOS} / M_{\star}^{BMA} \right)$'
res_abs_mag_label   = r'$M_r^{BMA}-M_r^{COSMOS}$'

def tree_variable_plot(x1,x2,x3,zlabel='$z_{cls}$'):
    cut = remove_nan(x2)
    idx = np.argsort(-1*x3)
    
    xmin, xmax = np.nanmin(np.hstack([x1,x2])), np.nanmax(np.hstack([x1,x2]))
    plt.plot([xmin,xmax],[xmin,xmax],'k--',lw=3)
    plt.scatter(x1[idx],x2[idx],c=x3[idx],cmap='RdBu',s=5,alpha=0.8)
    cbar = plt.colorbar()
    cbar.set_label(zlabel)
    
# put here below your plot functions

# Example of how to run
mass  = gal['mass']
mass_t= gal['masst']
magi  = gal['iobs']
magi_true = gal['Mi']
z = gal['redshift']
chisqr = gal['best_chisq']

mask = (z<=0.65)&(magi<23.5)&(mass_t>9.)
variables = [mass,mass_t,magi,magi_true,z]

# initiate class
b1 = bmaValidation()
b1.add_model('baseline',*variables)

# plot an identity plot
b1.plot_identity_mass_redshift(mask=mask)

# example: define new variables
b1.chisqr = np.log10(chisqr)
b1.plot_identity_mass_chisqr(mask=mask)
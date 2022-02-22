# Compute the cluster fraction for different labels
import numpy as np
import astropy.io.ascii as at

## stellar mass bin
msbins = np.arange(9.0,12.25,0.25)

## radii bins
rbins = np.arange(0.,3.25,0.25)

def compute_kde(x,weights,bw=0.1):
    Norm = np.nansum(weights)
    pdf = gaussian_kde(x, weights=weights, bw_method=bw)
    return pdf, Norm

def compute_fraction_kde(xvec,x,p1,p2,bw=None,eps=1e-9):
    pdf1, N1 = compute_kde(x, p1, bw=bw)
    pdf2, N2 = compute_kde(x, p1*p2, bw=bw)
    frac = N2*pdf2(xvec)/(N1*pdf1(xvec))
    return frac

def compute_fraction(prob,eps=1e-6):
    N1, N2 = np.nansum(prob), prob[~np.isnan(prob)].size
    frac = N1/N2
    frac_err = frac*np.sqrt(1/(N1+eps)+1/(N2+eps))
    return np.array([frac,frac_err])

def compute_fraction1(prob1, prob2, eps=1e-6):
    N1, N2 = np.nansum(prob1*prob2),prob1[~np.isnan(prob1)].size
    N3, N4 = np.nansum(prob1), np.nansum(prob2)
    
    frac = N1/N3
    #frac2= N1/N3
    #frac3= N1/N4
    
    frac_err = frac*np.sqrt(1/(N1+eps)+1/(N2+eps))
    frac_err = np.where(frac_err>1.,1.,frac_err)
    return np.array([frac,frac_err])

def compute_fraction2(prob1, prob2, eps=1e-6):
    N1, N2 = np.nansum(prob1*prob2), np.nansum(prob1)
    #N1, N2 = np.nansum(prob1*prob2),prob[~np.isnan(prob)].size
    
    frac = N1/N2
    frac_err = frac*np.sqrt(1/(N1+eps)+1/(N2+eps))
    return np.array([frac,frac_err])

def quenching_fraction_excess(fq1,fq2):
    dfrac = fq2[0]-fq1[0]
    qfe = dfrac/(1-fq1[0])
    qfe_err = np.abs(qfe)*np.sqrt((np.sqrt(fq2[1]**2+fq1[1]**2)/dfrac)**2 + (fq1[1]/(1-fq1[0]))**2)
    qfe_err = np.where(np.abs(qfe_err)>1.0,1.0,qfe_err)
    return np.array([qfe,qfe_err])

def make_bins(x,xbins):
    indices = []
    xmd     = 0.5*(xbins[1:]+xbins[:-1])
    for xl,xh in zip(xbins[:-1],xbins[1:]):
        w, = np.where((x<=xh)&(x>xl))
        indices.append(w)
    return indices,xmd

def chunks(ids1, ids2):
    """Yield successive n-sized chunks from data"""
    for id in ids2:
        w, = np.where( ids1==id )
        yield w
        
def save_output_matrix(var,outfile):
    nout = len(var)
    nsize = var[0].size
    out = np.zeros((2*nout, nsize) ,dtype=float)
    
    out[0] = var[0]
    for i,x in enumerate(var[1:]):
        ii = 2*i+1
        out[ii] = x[0]
        out[ii+1] = x[1]
    np.savetxt(outfile, out.T, fmt='%.5f')
    
def check_non_valid_number(x):
    w,  = np.where(x < 0)
    x[w] = np.nan
    return np.array(x)

label1 = ['Cluster', 'Infall', 'Interlopers']
label2 = ['Cluster+Infall', 'Cluster+Interlopers', 'Infall+Interlopers']

def plot_fraction_pannel(label, type):
    x = np.loadtx(file_base.format(label, type))
    
    fig = plt.figure(figsize=(12,4))

    plt.subplot(1, 2, 1)
    for i in range(3):
        ii = 2*i + 1
        plt.errorbar(x[0], x[ii], yerr=x[ii+1], label=label1[i], fmt='o')
    
    plt.legend(fontsize=12,loc=2)
    plt.xlabel(r'Log($M_\star/M_{\odot}$)',fontsize=16)
    plt.title('Quenching Fraction',fontsize=16)
    plt.ylim(-0.02,1.05)

    plt.subplot(1, 2, 2)
    for i in range(3, 6):
        ii = 2*i + 1
        plt.errorbar(x[0], x[ii], yerr=x[ii+1], label=label2[i], fmt='o')
    plt.legend(fontsize=12,loc=2)
    plt.xlabel(r'Log($M_\star/M_{\odot}$)',fontsize=18)
    plt.title('Quenching Fraction Excess',fontsize=16)

    plt.ylim(-0.02,1.05)
    
import numpy as np
import matplotlib.pyplot as plt


# class from here: http://nbviewer.ipython.org/gist/tillahoffmann/f844bce2ec264c1c8cb5
class gaussian_kde(object):
    """Representation of a kernel-density estimate using Gaussian kernels.
    Kernel density estimation is a way to estimate the probability density
    function (PDF) of a random variable in a non-parametric way.
    `gaussian_kde` works for both uni-variate and multi-variate data.   It
    includes automatic bandwidth determination.  The estimation works best for
    a unimodal distribution; bimodal or multi-modal distributions tend to be
    oversmoothed.
    Parameters
    ----------
    dataset : array_like
        Datapoints to estimate from. In case of univariate data this is a 1-D
        array, otherwise a 2-D array with shape (# of dims, # of data).
    bw_method : str, scalar or callable, optional
        The method used to calculate the estimator bandwidth.  This can be
        'scott', 'silverman', a scalar constant or a callable.  If a scalar,
        this will be used directly as `kde.factor`.  If a callable, it should
        take a `gaussian_kde` instance as only parameter and return a scalar.
        If None (default), 'scott' is used.  See Notes for more details.
    weights : array_like, shape (n, ), optional, default: None
        An array of weights, of the same shape as `x`.  Each value in `x`
        only contributes its associated weight towards the bin count
        (instead of 1).
    Attributes
    ----------
    dataset : ndarray
        The dataset with which `gaussian_kde` was initialized.
    d : int
        Number of dimensions.
    n : int
        Number of datapoints.
    neff : float
        Effective sample size using Kish's approximation.
    factor : float
        The bandwidth factor, obtained from `kde.covariance_factor`, with which
        the covariance matrix is multiplied.
    covariance : ndarray
        The covariance matrix of `dataset`, scaled by the calculated bandwidth
        (`kde.factor`).
    inv_cov : ndarray
        The inverse of `covariance`.
    Methods
    -------
    kde.evaluate(points) : ndarray
        Evaluate the estimated pdf on a provided set of points.
    kde(points) : ndarray
        Same as kde.evaluate(points)
    kde.pdf(points) : ndarray
        Alias for ``kde.evaluate(points)``.
    kde.set_bandwidth(bw_method='scott') : None
        Computes the bandwidth, i.e. the coefficient that multiplies the data
        covariance matrix to obtain the kernel covariance matrix.
        .. versionadded:: 0.11.0
    kde.covariance_factor : float
        Computes the coefficient (`kde.factor`) that multiplies the data
        covariance matrix to obtain the kernel covariance matrix.
        The default is `scotts_factor`.  A subclass can overwrite this method
        to provide a different method, or set it through a call to
        `kde.set_bandwidth`.
    Notes
    -----
    Bandwidth selection strongly influences the estimate obtained from the KDE
    (much more so than the actual shape of the kernel).  Bandwidth selection
    can be done by a "rule of thumb", by cross-validation, by "plug-in
    methods" or by other means; see [3]_, [4]_ for reviews.  `gaussian_kde`
    uses a rule of thumb, the default is Scott's Rule.
    Scott's Rule [1]_, implemented as `scotts_factor`, is::
        n**(-1./(d+4)),
    with ``n`` the number of data points and ``d`` the number of dimensions.
    Silverman's Rule [2]_, implemented as `silverman_factor`, is::
        (n * (d + 2) / 4.)**(-1. / (d + 4)).
    Good general descriptions of kernel density estimation can be found in [1]_
    and [2]_, the mathematics for this multi-dimensional implementation can be
    found in [1]_.
    References
    ----------
    .. [1] D.W. Scott, "Multivariate Density Estimation: Theory, Practice, and
           Visualization", John Wiley & Sons, New York, Chicester, 1992.
    .. [2] B.W. Silverman, "Density Estimation for Statistics and Data
           Analysis", Vol. 26, Monographs on Statistics and Applied Probability,
           Chapman and Hall, London, 1986.
    .. [3] B.A. Turlach, "Bandwidth Selection in Kernel Density Estimation: A
           Review", CORE and Institut de Statistique, Vol. 19, pp. 1-33, 1993.
    .. [4] D.M. Bashtannyk and R.J. Hyndman, "Bandwidth selection for kernel
           conditional density estimation", Computational Statistics & Data
           Analysis, Vol. 36, pp. 279-298, 2001.
    Examples
    --------
    Generate some random two-dimensional data:
    >>> from scipy import stats
    >>> def measure(n):
    >>>     "Measurement model, return two coupled measurements."
    >>>     m1 = np.random.normal(size=n)
    >>>     m2 = np.random.normal(scale=0.5, size=n)
    >>>     return m1+m2, m1-m2
    >>> m1, m2 = measure(2000)
    >>> xmin = m1.min()
    >>> xmax = m1.max()
    >>> ymin = m2.min()
    >>> ymax = m2.max()
    Perform a kernel density estimate on the data:
    >>> X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    >>> positions = np.vstack([X.ravel(), Y.ravel()])
    >>> values = np.vstack([m1, m2])
    >>> kernel = stats.gaussian_kde(values)
    >>> Z = np.reshape(kernel(positions).T, X.shape)
    Plot the results:
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111)
    >>> ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
    ...           extent=[xmin, xmax, ymin, ymax])
    >>> ax.plot(m1, m2, 'k.', markersize=2)
    >>> ax.set_xlim([xmin, xmax])
    >>> ax.set_ylim([ymin, ymax])
    >>> plt.show()
    """
    def __init__(self, dataset, bw_method=None, weights=None):
        self.dataset = np.atleast_2d(dataset)
        if not self.dataset.size > 1:
            raise ValueError("`dataset` input should have multiple elements.")
        self.d, self.n = self.dataset.shape

        if weights is not None:
            self.weights = weights / np.sum(weights)
        else:
            self.weights = np.ones(self.n) / self.n

        # Compute the effective sample size
        # http://surveyanalysis.org/wiki/Design_Effects_and_Effective_Sample_Size#Kish.27s_approximate_formula_for_computing_effective_sample_size
        self.neff = 1.0 / np.sum(self.weights ** 2)

        self.set_bandwidth(bw_method=bw_method)

    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.
        Parameters
        ----------
        points : (# of dimensions, # of points)-array
            Alternatively, a (# of dimensions,) vector can be passed in and
            treated as a single point.
        Returns
        -------
        values : (# of points,)-array
            The values at each point.
        Raises
        ------
        ValueError : if the dimensionality of the input points is different than
                     the dimensionality of the KDE.
        """
        from scipy.spatial.distance import cdist

        points = np.atleast_2d(points)

        d, m = points.shape
        if d != self.d:
            if d == 1 and m == self.d:
                # points was passed in as a row vector
                points = np.reshape(points, (self.d, 1))
                m = 1
            else:
                msg = "points have dimension %s, dataset has dimension %s" % (d, self.d)
                raise ValueError(msg)

        # compute the normalised residuals
        chi2 = cdist(points.T, self.dataset.T, 'mahalanobis', VI=self.inv_cov) ** 2
        # compute the pdf
        result = np.sum(np.exp(-.5 * chi2) * self.weights, axis=1) / self._norm_factor

        return result

    __call__ = evaluate

    def scotts_factor(self):
        return np.power(self.neff, -1. / (self.d + 4))

    def silverman_factor(self):
        return np.power(self.neff * (self.d + 2.0) / 4.0, -1. / (self.d + 4))

    #  Default method to calculate bandwidth, can be overwritten by subclass
    covariance_factor = scotts_factor

    def set_bandwidth(self, bw_method=None):
        """Compute the estimator bandwidth with given method.
        The new bandwidth calculated after a call to `set_bandwidth` is used
        for subsequent evaluations of the estimated density.
        Parameters
        ----------
        bw_method : str, scalar or callable, optional
            The method used to calculate the estimator bandwidth.  This can be
            'scott', 'silverman', a scalar constant or a callable.  If a
            scalar, this will be used directly as `kde.factor`.  If a callable,
            it should take a `gaussian_kde` instance as only parameter and
            return a scalar.  If None (default), nothing happens; the current
            `kde.covariance_factor` method is kept.
        Notes
        -----
        .. versionadded:: 0.11
        Examples
        --------
        >>> x1 = np.array([-7, -5, 1, 4, 5.])
        >>> kde = stats.gaussian_kde(x1)
        >>> xs = np.linspace(-10, 10, num=50)
        >>> y1 = kde(xs)
        >>> kde.set_bandwidth(bw_method='silverman')
        >>> y2 = kde(xs)
        >>> kde.set_bandwidth(bw_method=kde.factor / 3.)
        >>> y3 = kde(xs)
        >>> fig = plt.figure()
        >>> ax = fig.add_subplot(111)
        >>> ax.plot(x1, np.ones(x1.shape) / (4. * x1.size), 'bo',
        ...         label='Data points (rescaled)')
        >>> ax.plot(xs, y1, label='Scott (default)')
        >>> ax.plot(xs, y2, label='Silverman')
        >>> ax.plot(xs, y3, label='Const (1/3 * Silverman)')
        >>> ax.legend()
        >>> plt.show()
        """
        from six import string_types
        if bw_method is None:
            pass
        elif bw_method == 'scott':
            self.covariance_factor = self.scotts_factor
        elif bw_method == 'silverman':
            self.covariance_factor = self.silverman_factor
        elif np.isscalar(bw_method) and not isinstance(bw_method, string_types):
            self._bw_method = 'use constant'
            self.covariance_factor = lambda: bw_method
        elif callable(bw_method):
            self._bw_method = bw_method
            self.covariance_factor = lambda: self._bw_method(self)
        else:
            msg = "`bw_method` should be 'scott', 'silverman', a scalar " \
                  "or a callable."
            raise ValueError(msg)

        self._compute_covariance()

    def _compute_covariance(self):
        """Computes the covariance matrix for each Gaussian kernel using
        covariance_factor().
        """
        self.factor = self.covariance_factor()
        # Cache covariance and inverse covariance of the data
        if not hasattr(self, '_data_inv_cov'):
            # Compute the mean and residuals
            _mean = np.sum(self.weights * self.dataset, axis=1)
            _residual = (self.dataset - _mean[:, None])
            # Compute the biased covariance
            self._data_covariance = np.atleast_2d(np.dot(_residual * self.weights, _residual.T))
            # Correct for bias (http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_covariance)
            self._data_covariance /= (1 - np.sum(self.weights ** 2))
            self._data_inv_cov = np.linalg.inv(self._data_covariance)

        self.covariance = self._data_covariance * self.factor ** 2
        self.inv_cov = self._data_inv_cov / self.factor ** 2
        self._norm_factor = np.sqrt(np.linalg.det(2 * np.pi * self.covariance))  # * self.n
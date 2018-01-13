'''Functions in this module will deal with randomisation and statisitics on iCLIP 
data'''

import pandas

def count_ratio(x, test, control):
    '''Calculate the ratio of the sum to two columns of a DataFrame.
    
    Parameters
    ----------
    x : pandas.DataFrame
        DataFrame containing the columns to be used
    test : str
        Name of column to be used as the numerator in the ratio
    control : str
        Name of column to be used as the demoninator in the ratio

    Returns
    -------
    float
        Ratio of the sum of counts in column ``test`` divided by sum of
        counts in column ``control```. If the sum of the control column
        is 0, returns numpy.Inf
    '''
    
    try:
        return float(x[test].sum())/x[control].sum()
    except ZeroDivisionError:
        return pandas.np.inf
    
def bootstrap(x, bootstraps, fun, *args, **kwargs):
    '''Calculate bootstraps a function operating on a Series.

    Draw samples with replacement from `x` with replacement and apply `fun`
    . Do this `bootstrap` times.

    Parameters
    ----------
    x : pandas.Series or pandas.DataFrame
        `pandas` object to bootstrap and apply function to.
    bootstraps : int
        How many bootstrap samples to draw
    fun : func
        Function to apply to the bootstraps resamples of `x`. `x` will be
        passed as the first argument.
    *args, **kwargs : 
        Further arguments and keyword arguments passed on to `fun`.

    Returns
    -------
    pandas.Series of same type as ``fun(x)``
         Each entry in the Series will be the result of applying
         `fun` to a bootstrap sample from `x`.

    See Also
    --------
    boot_ci : Calculate confidence intervals by bootstrapping

    '''
    
    def _inner(y):
        xsample = x.sample(n=x.shape[0], replace=True)
        return fun(xsample, *args, **kwargs)
    return pandas.Series(range(bootstraps+1)).apply(_inner)

def boot_ci(x, quantiles=(0.05,0.95), *args,**kwargs):
    '''Calculate confidence intervals on the results of function
    using bootstrapping. 

    Parameters
    ----------
    x : pandas.Series or pandas.DataFrame
        Object containing data to bootstrap from
    quantiles : tuple of (float, float), optional
       quantiles to use in calculating confidence interval. For example,
       for a 95% confidence interval, supply (0.025, 0.975)
    *args, **kwargs :
       further arguments and keyword arguements passed to :func:`bootstrap`
       (and likely onwards to the function called by bootstrap).

    Returns
    -------
    pandas.Series
        Series with the specified quantiles from the bootstraps

    
    Other Parameters
    ----------------
    fun : func
        function to be applied to bootstrap samples of x
    bootstraps : int
        number of bootstrap samples to use_args

    See Also
    --------
    bootstrap : draw bootstrap samples and apply a function to eaach
    ratio_and_ci : Uses this function to calculate confidence intervals on
                   the ratio of two columns. 
    
    
    '''

    return bootstrap(x,*args, **kwargs).quantile(quantiles)


def ratio_and_ci(test, control, quantiles=(0.05, 0.95), bootstraps=1000):
    '''Calculating ratio between the sum of two columns and add a
    confidence interval calculated using bootstrapping.

    Parameters
    ----------
    test : sequence or pandas.Series
        set of numbers to use in calculating the numerator of the ratio.
    control : sequence or pandas.Series
        set of numbers to use in calculating the domoninator of the ratio.
    quantiles : tuple of (float, float)
        quantiles to use in calculating confidence interval. For example,
        for a 95% confidence interval, supply (0.025, 0.975).
    bootstraps : int, optional
        Number of bootstrap samples to draw.

    Returns
    -------
    pandas.Series
        Series with three entries - ratio, and two quantiles, formed by 
        appending q to the request quantiles: e.g. for default quantiles
        the index will be ["ratio", "q0.05", "q0.95"]

    See also
    --------
    bootstrap, boot_ci

    Notes
    -----
    Ideally, the test and control will be of the same length. When bootstraps
    are drawn, they are drawn in pairs, one from test, one from control.
    If test and control are similarly index series (e.g. counts from the same
    transcripts or exons) then this pairing is maintained. This hopefully 
    helps to counteract some of the effects of correlation between the 
    data sets.


    '''
    
    x = pandas.DataFrame({"test": test, "control": control})
    x_ci = boot_ci(x, quantiles=quantiles,
                   bootstraps=bootstraps,
                   fun=count_ratio,
                   test="test",
                   control="control")
    ratio = count_ratio(x, "test", "control")
    try:
        results={"q%s" % q: x for q, x in x_ci.iteritems()}
    except AttributeError:
        results = {"q%s" % quantiles: x_ci}

    results["ratio"] = ratio
    
    return pandas.Series(results)



'''Functions in this module will deal with randomisation and statisitics on iCLIP 
data'''

import pandas

def count_ratio(x, test, control):
    '''Sum counts in the test_FLAG_union column of x and divide
    by the sum of the control column of x. If the sum of the control
    column is 0, return Inf'''
    try:
        return float(x[test].sum())/x[control].sum()
    except ZeroDivisionError:
        return pandas.np.inf
    
def bootstrap(x, bootstraps, fun, *args, **kwargs):
    '''Calculate bootstrap bootstraps of the function fun
    from x by sampling with replacement from the rows of
    x and computing `fun` on the results. args and kwargs 
    passed on to fun'''
    def _inner(y):
        xsample = x.sample(n=x.shape[0], replace=True)
        return fun(xsample, *args, **kwargs)
    return pandas.Series(range(bootstraps+1)).apply(_inner)

def boot_ci(x, quantiles=(0.05,0.95), *args,**kwargs):
    '''Calculate the quantiles for bootstraps of x
    args and kwargs are passed on to bootstraps'''
    return bootstrap(x,*args, **kwargs).quantile(quantiles)


def ratio_and_ci(test, control, quantiles=(0.05, 0.95), bootstraps=1000):
    '''Automates calculating ratio between a test and a control for each
    set in groupby and add a 90% confidence interval to it. Groupby
    can be anything that would work with x.groupby(groupby). Uses 1000
    bootstraps
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



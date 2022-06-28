import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm


def func(x, xs, log_ys):
    mu,sigma,offset = x
    return np.sum(
        [(np.log(norm.pdf(x,loc=mu,scale=sigma))+np.log(offset)-log_y)**2 \
            for x,log_y in zip(xs,log_ys)])

def fit_normal(vals, init, bounds, sgn):
    vals = np.sort(vals)
    vals = vals[:len(vals)//2] if method==1 else vals[len(vals)//2:]
    xs, ys = np.unique(vals, return_counts=True)
    log_ys = np.log(ys/sum(ys))
    mu,sigma,offset = minimize(func, x0=init, args=(xs,log_ys), bounds=bounds).x
    print('mu:', mu, 'stdv:', sigma, 'offset:', offset, 'thresh:',mu+sgn*2*sigma)
    return mu,sigma

def get_threshold(sketches, seqs, est_overlap, **args):
    vals = []
    n = len(seqs)
    for i in np.random.choice(np.arange(n), size=n//10, replace=False):
        for j in np.random.choice(np.arange(n), size=n//20, replace=False):
            if i==j: continue
            _, arr = est_overlap(sketches[i], sketches[j], seqs[i], seqs[j])
            vals.extend(arr)
    mu,sigma = fit_normal(vals, **args)
    thresh = mu + args['sgn']*2*sigma
    return thresh
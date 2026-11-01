# Poisson GLM by iteratively reweighted least squares in numpy, replacing statsmodels' GLM.fit for the many bootstrap fits of
# the fire-activity regression (same model, same maximum-likelihood estimate; compare with statsmodels by running cloudfree_analysis.py check with FASTGLM=0 and 1).
import numpy as np
def poisson_irls(y, X, off, tol=1e-8, maxit=100):
    """Same stopping rule as statsmodels' IRLS: the deviance changes by less than tol (absolute) between iterations."""
    y = np.asarray(y, float); X = np.asarray(X, float); off = np.zeros(len(y)) if off is None else np.asarray(off, float)
    mu = (y + y.mean()) / 2; eta = np.log(mu); dev = np.inf
    ylogy = np.where(y > 0, y * np.log(np.where(y > 0, y, 1)), 0.0)
    for it in range(maxit):
        z = eta - off + (y - mu) / mu
        Xw = X * mu[:, None]; H = X.T @ Xw
        if np.linalg.cond(H) < 1e10: b = np.linalg.solve(H, Xw.T @ z)          # weighted normal equations when well conditioned
        else: w = np.sqrt(mu); b = np.linalg.lstsq(X * w[:, None], z * w, rcond=None)[0]   # otherwise minimum-norm least squares
        eta = X @ b + off; mu = np.exp(eta)
        dn = 2 * np.sum(ylogy - y * eta - (y - mu))
        if abs(dn - dev) < tol: return b
        dev = dn
    raise RuntimeError('IRLS did not converge')
class _Res:
    def __init__(self, params): self.params = params
def install(sm):
    import statsmodels.genmod.families as fam
    class GLM:
        def __init__(self, endog, exog, family=None, offset=None):
            assert isinstance(family, fam.Poisson); self.args = (endog, exog, offset)
        def fit(self): return _Res(poisson_irls(*self.args))
    sm.GLM = GLM

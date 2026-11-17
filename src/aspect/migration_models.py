"""
Model comparison on the migration front (revision version).

Migration models (front grows away from the origin):
  ramp       D = v * min(T, tau)            linear migration that stops at tau   (k = 2)
  sqrt_ramp  D = a * sqrt(min(T, tau))      diffusion-like migration that stops  (k = 2)
  linear     D = v * T                      linear migration, no stop            (k = 1)
  sqrt       D = a * sqrt(T)                diffusion, no stop                   (k = 1)
  (tau must leave >= min_event_frac of the events inside the migration segment)

Non-migration models:
  const      D = c                          no growth                            (k = 1)
  step       D = c1 | c2                    sudden jump of the front             (k = 3)
  log        D = a + b ln T, b > 0          aftershock-like expansion            (k = 2)

IC = n_eff * ln(SSE / n) + k * ln(n_eff), with n_eff = n * step / window
(overlapping windows are not independent).
"""
import numpy as np

MIGRATION_MODELS = ('ramp', 'sqrt_ramp', 'linear', 'sqrt')
SHAPE = {'ramp': lambda x: x, 'linear': lambda x: x, 'sqrt_ramp': np.sqrt, 'sqrt': np.sqrt}


def fit_models(T, D, w, step, n_events, min_event_frac=0.3, pure_models=True, log_positive=True):
    """Fit all models; returns a dict model -> {ic, parameters}, or None if the front is too short."""
    ok = (T > 0) & (D > 0)
    T, D = T[ok], D[ok]
    n = len(T)
    # smallest index the stop time may take so that the segment holds >= min_event_frac of the events
    min_pts = max(4, int(np.ceil((min_event_frac * n_events - w) / step)) + 1)
    if n < 8 or min_pts > n:
        return None
    n_eff = max(n * step / w, 3.0)
    ic = lambda sse, k: n_eff * np.log(max(sse, 1e-12 * n) / n) + k * np.log(n_eff)

    M = {'const': dict(ic=ic(np.sum((D - D.mean()) ** 2), 1))}

    # aftershock-like expansion D = a + b ln T
    X = np.c_[np.ones(n), np.log(T)]
    b = np.linalg.lstsq(X, D, rcond=None)[0]
    good = b[1] > 0 or not log_positive
    M['log'] = dict(ic=ic(np.sum((D - X @ b) ** 2), 2) if good else np.inf, a=b[0], b=b[1])

    # jump: best split of the front into two constant levels
    cs, cs2 = np.cumsum(D), np.cumsum(D * D)
    best = None
    for i in range(2, n - 1):
        sse = cs2[i - 1] - cs[i - 1] ** 2 / i + (cs2[-1] - cs2[i - 1]) - (cs[-1] - cs[i - 1]) ** 2 / (n - i)
        if best is None or sse < best[0]:
            best = (sse, i)
    M['step'] = dict(ic=ic(best[0], 3), split=T[best[1]], c1=D[:best[1]].mean(), c2=D[best[1]:].mean())

    # migrations that stop: grid search over the stop time
    js = np.arange(min_pts - 1, n)
    for name in ('ramp', 'sqrt_ramp'):
        M[name] = dict(ic=np.inf)
        Xr = SHAPE[name](np.minimum(T[None, :], T[js][:, None]))
        a = (Xr @ D) / np.einsum('cn,cn->c', Xr, Xr)
        sse = np.where(a > 0, ((D - a[:, None] * Xr) ** 2).sum(1), np.inf)
        j = int(np.argmin(sse))
        if np.isfinite(sse[j]):
            M[name] = dict(ic=ic(sse[j], 2), rate=float(a[j]), tau=float(T[js[j]]))

    # migrations without a stop (one parameter)
    for name in ('linear', 'sqrt'):
        M[name] = dict(ic=np.inf)
        if pure_models:
            x = SHAPE[name](T)
            a = float(x @ D / (x @ x))
            if a > 0:
                M[name] = dict(ic=ic(np.sum((D - a * x) ** 2), 1), rate=a, tau=float(T[-1]))
    M['_T'], M['_D'] = T, D
    return M


def predict(name, m, T):
    """Model curve on times T."""
    if name in MIGRATION_MODELS:
        return m['rate'] * SHAPE[name](np.minimum(T, m['tau']))
    if name == 'const':
        return np.full_like(T, np.nan, dtype=float)
    if name == 'log':
        return m['a'] + m['b'] * np.log(T)
    if name == 'step':
        return np.where(T < m['split'], m['c1'], m['c2'])


def migration_phase_r2(M, name):
    """Score the selected growth fit only at front centres 0 < t <= tau.

    Use the same unsmoothed distance front and fitted amplitude as model
    comparison, including for diffusion. Never optimize a second prefix,
    refit the amplitude, transform distance to distance squared, or clip R2.
    """
    if name is None:
        return dict(r2=-np.inf, r2_n_points=0, r2_tau=np.nan)
    model = M[name]
    T, D = M['_T'], M['_D']
    mask = (T > 0) & (T <= model['tau']) & np.isfinite(D)
    y = D[mask]
    score = -np.inf
    if len(y) >= 4:
        sst = float(np.sum((y - y.mean()) ** 2))
        if sst > 1e-12:
            sse = float(np.sum((y - predict(name, model, T[mask])) ** 2))
            score = 1.0 - sse / sst
    return dict(r2=score, r2_n_points=int(mask.sum()), r2_tau=float(model['tau']))


def compare(M, step_margin=3.0, log_margin=6.0):
    """
    dIC = min(IC_const, IC_step + step_margin, IC_log + log_margin) - IC_best_migration.
    The jump / log expansion vetoes a migration only if it fits better by more than its margin.
    """
    mig = min(MIGRATION_MODELS, key=lambda k: M[k]['ic'])
    terms = {'const': M['const']['ic'], 'step': M['step']['ic'] + step_margin, 'log': M['log']['ic'] + log_margin}
    veto = min(terms, key=terms.get)
    return dict(mig_best=mig if np.isfinite(M[mig]['ic']) else None, veto=veto,
                dic=float(terms[veto] - M[mig]['ic']),
                jump_adv=float(M[mig]['ic'] - M['step']['ic']),   # > 0: jump fits better
                log_adv=float(M[mig]['ic'] - M['log']['ic']))     # > 0: log t fits better

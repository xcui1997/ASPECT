"""
Migration test (revision version): is there a migration, and how well are its parameters constrained?

For one sequence:
  1. front + model comparison (step 1) and migration-phase R2 (step 2) on the selected growth segment ("main run");
  2. the same on 18 perturbed versions (random 90 % subsamples, random subsets retaining 90 % of events);
  3. vote = share of the 19 runs passing both steps; migration if vote >= vote_threshold;
  4. velocity stability uses jointly passing fits from these same runs; no new bootstrap by default.
All parameters come from config['migration_test_params'].
"""
import numpy as np
import matplotlib.pyplot as plt

from .migration_front import trim_early_events, build_front, distance_to_origin
from .migration_models import fit_models, predict, compare, SHAPE, migration_phase_r2

MODEL_NAME = {'ramp': 'linear, then stop', 'sqrt_ramp': 'diffusion, then stop', 'linear': 'linear', 'sqrt': 'diffusion',
              'const': 'stationary', 'step': 'jump', 'log': 'log t expansion'}
REASON = {'step': 'N1 jump fits better', 'log': 'N2 log t expansion fits better', 'const': 'N3 stationary fits better'}


class MigrationTester:

    def __init__(self, config):
        p = config['migration_test_params']
        self.params = p
        # These fields state the confirmed scientific representation explicitly.
        required = {'front_definition':'window_center_interpolated_distance', 'percentile_method':'linear',
                    'growth_coverage':'window_index_approximation', 'front_deduplicate_events':False}
        for key, value in required.items():
            if p.get(key, value) != value:
                raise ValueError(f'{key} must be {value!r} for the current workflow')
        self.pct = p.get('distance_percentile', 70)
        self.n_first = p.get('origin_n_first', 10)
        self.window_ratio = p.get('window_ratio', 0.05)
        self.step_ratio = p.get('step_ratio', 0.4)
        self.min_event_frac = p.get('min_event_frac', 0.3)
        self.pure_models = p.get('pure_models', True)
        self.log_positive = p.get('log_positive', True)
        self.step_margin = p.get('step_margin', 3.0)
        self.log_margin = p.get('log_margin', 6.0)
        self.dic_threshold = p.get('dic_threshold', 0.0)
        self.r2 = p.get('r2', {})
        self.r2_threshold = self.r2.get('threshold', 0.5)
        self.pert = p.get('perturbation', {})
        if self.pert.get('n_subsample', 18) != 18:
            raise ValueError('Current migration protocol requires exactly 18 random 90% subsets plus the baseline.')
        if self.pert.get('origin_n_first', []) or self.pert.get('percentiles', []):
            raise ValueError('Alternative origins and percentiles are separate sensitivity tests, not members of the 19-run vote.')
        if not np.isclose(self.pert.get('drop_rate', 0.1), 0.1):
            raise ValueError('Current migration protocol retains 90% of events in every subset.')
        self.vote_threshold = p.get('vote_threshold', 0.5)
        self.reason_vote_floor = p.get('reason_vote_floor', 0.3)
        self.early = p.get('early_screen', {})
        self.unc = p.get('uncertainty', {})

    # ------------------------------------------------------------------ one run
    def _one_run(self, coords, time, pct, n_first):
        """Compare models, then score the selected fit on its migration phase only."""
        T, D, w, step = build_front(coords, time, pct, n_first, self.window_ratio, self.step_ratio)
        M = fit_models(T, D, w, step, len(time), self.min_event_frac, self.pure_models, self.log_positive)
        if M is None:
            return dict(valid=False, r2=-np.inf, r2_n_points=0, r2_tau=np.nan,
                        dic=-np.inf, passed=False, mig_best=None, models=None)
        out = compare(M, self.step_margin, self.log_margin)
        out.update(valid=True, models=M, **migration_phase_r2(M, out['mig_best']))
        out['passed'] = out['dic'] > self.dic_threshold and out['r2'] > self.r2_threshold
        return out

    def _perturbations(self, coords, time, rng):
        """Exactly 18 random 90 % subsamples; origin rule and front percentile stay fixed."""
        n = len(time)
        keep = int((1 - self.pert.get('drop_rate', 0.1)) * n)
        versions = []
        for _ in range(self.pert.get('n_subsample', 18)):
            idx = np.sort(rng.choice(n, keep, replace=False))
            versions.append((coords[idx], time[idx], self.pct, self.n_first))
        return versions

    # ------------------------------------------------------------------ full test
    def analyze(self, coords, time, mag=None, seed=0, per_run_initial_trim_fraction=None):
        """Migration test for one sequence (coords in km, time in days). Returns a result dict or None."""
        if self.early.get('enabled', True):
            coords, time, mag, start = trim_early_events(coords, time, mag, self.early.get('lower_percentile', 10),
                                                         self.early.get('upper_percentile', 90), self.early.get('multiplier', 1.0))
        else:
            coords, time, mag, start = trim_early_events(coords, time, mag, min_keep=10 ** 9)   # sort only
        if len(time) < self.params.get('min_events_required', 50):
            return None
        rng = np.random.default_rng(seed)

        # Optional historical case-specific sensitivity procedure: draw each
        # subset from the full screened sequence, THEN remove its earliest
        # floor(fraction * n) events. This is separate from early_screen.
        # The automatic production default is zero (no case-specific trim).
        fraction = self.params.get('per_run_initial_trim_fraction', 0.0) if per_run_initial_trim_fraction is None else per_run_initial_trim_fraction
        if not 0 <= fraction < 1:
            raise ValueError('per_run_initial_trim_fraction must be in [0, 1)')
        versions = [(coords, time, self.pct, self.n_first)] + self._perturbations(coords, time, rng)
        trimmed_versions = []
        for c, t, pct, n_first in versions:
            cut = int(np.floor(fraction * len(t)))
            c, t = c[cut:], t[cut:]
            if cut and len(t):
                t = t - t[0]
            trimmed_versions.append((c, t, pct, n_first))
        # The event-count gate also applies after any optional full-data trim;
        # random subsets themselves need not contain 50 events.
        if len(trimmed_versions[0][1]) < self.params.get('min_events_required', 50):
            return None
        main = self._one_run(*trimmed_versions[0])
        runs = [main] + [self._one_run(*v) for v in trimmed_versions[1:]]
        cut_main = int(np.floor(fraction * len(time)))
        if cut_main:
            coords, time = trimmed_versions[0][:2]
            mag = None if mag is None else mag[cut_main:]
            start += cut_main
        # Retain the selected growth model for each run; summarize speed only over joint passes.
        run_parameters = []
        for version, run in enumerate(runs):
            model = run['mig_best']
            fit = run['models'][model] if run['valid'] and model else None
            tau = fit['tau'] if fit else np.nan
            distance = float(predict(model, fit, np.array([tau]))[0]) if fit else np.nan
            run_parameters.append(dict(version=version, valid=run['valid'], joint_pass=run['passed'],
                r2=run['r2'], dic=run['dic'], migration_model=model, duration=tau,
                distance=distance, velocity=distance/tau if tau > 0 else np.nan))
        vote = float(np.mean([r['passed'] for r in runs]))
        is_mig = vote >= self.vote_threshold

        # reason for 'no migration' (checked in this order)
        if is_mig:
            reason = 'migration'
        elif not main['valid']:
            reason = 'N0 too short to fit'
        elif vote >= self.reason_vote_floor or main['passed']:
            reason = f'N5 fewer than {int(np.ceil(self.vote_threshold * len(runs)))} of {len(runs)} runs pass'
        elif main['dic'] > self.dic_threshold:
            reason = f'N4 migration fits better but R2 <= {self.r2_threshold:g}'
        else:
            reason = REASON[main['veto']]

        res = dict(coords=coords, time=time, mag=mag, start_index=start, n_events=len(time),
                   migration=is_mig, reason=reason, vote=vote, main=main)
        if main['mig_best'] is not None:
            m = main['models'][main['mig_best']]
            dist = float(predict(main['mig_best'], m, np.array([m['tau']]))[0])
            res.update(duration=m['tau'], distance=dist, velocity=dist / m['tau'] if m['tau'] > 0 else np.nan)
        if is_mig and self.unc.get('enabled', False):
            res['uncertainty'] = self.fixed_window_uncertainty(coords, time, main['mig_best'], res['duration'],
                                                               np.random.default_rng(seed + self.unc.get('seed_offset', 777)))
        res['runs'] = run_parameters
        return res

    def fixed_window_uncertainty(self, coords, time, model, tau, rng):
        """
        Keep the migration model and its stop time tau; refit only the amplitude on n_draws random
        subsamples (drop_rate of the events removed). distance = model(tau), velocity = distance / tau.
        This is the sampling precision; it does not include full membership, model or location uncertainty.
        """
        n = len(time)
        keep = int(round((1 - self.unc.get('drop_rate', 0.1)) * n))
        dists = []
        for _ in range(self.unc.get('n_draws', 100)):
            idx = np.sort(rng.choice(n, keep, replace=False))
            T, D, _, _ = build_front(coords[idx], time[idx], self.pct, self.n_first, self.window_ratio, self.step_ratio)
            ok = (T > 0) & (D > 0)
            if ok.sum() < 5:
                continue
            x = SHAPE[model](np.minimum(T[ok], tau))
            dists.append(float(x @ D[ok] / (x @ x)) * float(SHAPE[model](tau)))
        q = np.percentile(dists, [2.5, 16, 50, 84, 97.5])
        out = {f'distance_p{k}': v for k, v in zip(('025', '16', '50', '84', '975'), q)}
        out.update({f'velocity_p{k}': v / tau for k, v in zip(('025', '16', '50', '84', '975'), q)})
        out['rel_halfwidth_68'] = (q[3] - q[1]) / 2 / q[2]
        out['rel_halfwidth_95'] = (q[4] - q[0]) / 2 / q[2]
        return out

    # ------------------------------------------------------------------ summary and figure
    def summary_row(self, res):
        """Flat dict for the output table."""
        main = res['main']
        row = dict(n_events=res['n_events'], start_index=res['start_index'], migration=res['migration'],
                   reason=res['reason'], vote=res['vote'], dic=main['dic'], r2=main['r2'],
                   r2_n_points=main['r2_n_points'], r2_tau=main['r2_tau'],
                   r2_definition='selected_model_window_center_phase_distance',
                   migration_model=main.get('mig_best'), non_migration_model=main.get('veto'),
                   jump_adv=main.get('jump_adv', np.nan), log_adv=main.get('log_adv', np.nan),
                   duration=res.get('duration', np.nan), distance=res.get('distance', np.nan),
                   velocity=res.get('velocity', np.nan))
        row.update(res.get('uncertainty', {}))
        # Use sample standard deviation (ddof=1) over joint passes, not a confidence interval.
        velocities = np.array([r['velocity'] for r in res['runs'] if r['joint_pass']])
        velocities = velocities[np.isfinite(velocities)]
        mean = float(velocities.mean()) if len(velocities) else np.nan
        sd = float(velocities.std(ddof=1)) if len(velocities)>1 else np.nan
        row.update(velocity_n=len(velocities), velocity_mean=mean, velocity_sd=sd,
                   velocity_variance=sd**2, velocity_cv=sd/mean if mean>0 else np.nan,
                   passing=sum(r['joint_pass'] for r in res['runs']))
        return row

    def plot_report(self, res, title='', save_path=None):
        """Left: distance vs time with front and best models; right: map on the two main axes coloured by time."""
        coords, time, main = res['coords'], res['time'], res['main']
        dist, origin = distance_to_origin(coords, self.n_first)
        T, D, _, _ = build_front(coords, time, self.pct, self.n_first, self.window_ratio, self.step_ratio)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), gridspec_kw=dict(width_ratios=[1.5, 1]))
        ax1.scatter(time, dist, s=4, c='#c9c8c2', edgecolors='none', label='events')
        ax1.plot(T, D, color='#3a3936', lw=1, label=f'P{self.pct} front')
        M = main.get('models')
        if M is not None:
            tt = np.linspace(M['_T'].min(), M['_T'].max(), 300)
            if main['mig_best']:
                ax1.plot(tt, predict(main['mig_best'], M[main['mig_best']], tt), color='#2a78d6', lw=1.8,
                         label=f"migration: {MODEL_NAME[main['mig_best']]}")
            nb = main['veto']
            yy = np.full_like(tt, M['_D'].mean()) if nb == 'const' else predict(nb, M[nb], tt)
            ax1.plot(tt, yy, color='#eb6834', lw=1.4, ls='--', label=f'non-migration: {MODEL_NAME[nb]}')
        ax1.set_xlabel('Time since first event (days)'); ax1.set_ylabel('3D distance to origin (km)')
        ax1.set_ylim(0, max(np.percentile(dist, 98), np.max(D)) * 1.1); ax1.legend(fontsize=7, frameon=False)
        label = 'migration' if res['migration'] else res['reason']
        ax1.set_title(f"{title}\n{label} | vote={res['vote']:.2f}  dIC={main['dic']:.1f}  R2={main['r2']:.2f}", fontsize=9)

        c = coords - origin
        _, _, vt = np.linalg.svd(c - c.mean(0), full_matrices=False)
        p = c @ vt[:2].T
        sc = ax2.scatter(p[:, 0], p[:, 1], s=5, c=time, cmap='viridis', edgecolors='none')
        ax2.plot(0, 0, '*', color='#d03b3b', ms=10)
        ax2.set_aspect('equal', 'datalim'); ax2.set_xlabel('Main axis 1 (km)'); ax2.set_ylabel('Main axis 2 (km)')
        fig.colorbar(sc, ax=ax2, label='Time (days)')
        fig.tight_layout()
        if save_path is not None:
            fig.savefig(save_path, dpi=150); plt.close(fig)
        return fig

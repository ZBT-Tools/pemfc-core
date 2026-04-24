"""Coupling scheduler between the two-phase sub-model and the outer
current-density fixed-point iteration.

The two-phase (saturation / evaporation) dynamics and the outer current-
density update form a coupled fixed-point system. For some operating points
the coupling is well-behaved and the error decays monotonically; for others
(typically mid-voltage limiting-current regions) tight coupling produces a
sustained limit cycle whose error floor never reaches the convergence target.

Operator splitting — updating the two-phase sub-system only every N main
iterations — decouples the two time scales and breaks the limit cycle.
However, starting with tight coupling (N=1) and escalating after detecting a
plateau does NOT work: the early N=1 phase establishes a bad-basin fixed
point that the later N=5 phase cannot escape. The successful strategy is the
opposite — start pessimistic at N=n_high, and de-escalate to N=n_low once
the outer error has dropped below the empirical limit-cycle floor. Easy
points (monotonic decay) cross this floor and enjoy the faster N=n_low
iterations for the rest of the run; hard points (true limit cycle) stay
pinned above the floor and remain at N=n_high throughout.
"""


class TwoPhaseCouplingScheduler:
    """Picks the two-phase update frequency (N) based on outer-error history.

    Starts at ``n_high`` and de-escalates once to ``n_low`` when the outer
    error drops below ``error_threshold`` — a level empirically chosen to
    sit just below the N=n_high limit-cycle floor, so only genuinely
    converging operating points trigger de-escalation. De-escalation is
    one-shot per operating point; call :meth:`reset` at the start of each
    new operating point.
    """

    def __init__(self, n_high: int = 5, n_low: int = 1,
                 error_threshold: float = 1e-5):
        self.n_high = n_high
        self.n_low = n_low
        self.error_threshold = error_threshold
        self.reset()

    def reset(self):
        """Called at the start of each new operating point."""
        self._n = self.n_high

    @property
    def n(self) -> int:
        return self._n

    def should_update_two_phase(self, iteration: int) -> bool:
        """True when the two-phase sub-system should be updated this step."""
        return iteration % self._n == 0

    def register_error(self, error: float):
        """Report the latest outer-iteration error. May de-escalate N."""
        if self._n == self.n_low:
            return
        # Ignore non-positive errors — happens on the first iteration when
        # current_density and current_density_old are both at their initial
        # values, giving a spurious RRMSE of 0. A real post-update error is
        # strictly positive.
        if error <= 0.0:
            return
        if error < self.error_threshold:
            self._n = self.n_low

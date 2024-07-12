from dataclasses import dataclass


@dataclass
class GlobalState:
    iteration: int
    error: float

    def __init__(self, iteration=0, error=1e8):
        self.iteration = iteration
        self.error = error


global_state = GlobalState()


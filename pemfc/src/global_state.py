from dataclasses import dataclass


@dataclass
class GlobalState:
    iteration: int
    error: float
    base_directory: str
    max_iteration: int


global_state = GlobalState(
    iteration=0,
    error=1e8,
    base_directory='.',
    max_iteration=1000)

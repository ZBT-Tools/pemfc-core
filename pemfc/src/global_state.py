from dataclasses import dataclass


@dataclass
class GlobalState:
    iteration: int
    error: float
    base_directory: str


global_state = GlobalState(
    iteration=0,
    error=1e8,
    base_directory='.')

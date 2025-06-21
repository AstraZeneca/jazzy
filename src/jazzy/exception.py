"""Anything about exceptions and their handling."""
# src/jazzy/exception.py
from functools import wraps


class JazzyError(Exception):
    """Raised when Jazzy cannot calculate results."""

    pass


class KallistoError(Exception):
    """Raised when something goes wrong inside kallisto."""

    pass


class NegativeLonePairsError(Exception):
    """Raised when a molecule contains negative lone pairs."""

    pass


class EmbeddingError(Exception):
    """Raised when the embedding fails for a molecule."""

    pass


def exception_handling(func):
    """Catches core exceptions and raises a JazzyError."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except (NegativeLonePairsError, KallistoError) as e:
            raise JazzyError(str(e))

    return wrapper

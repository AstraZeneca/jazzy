"""Data structure helpers. Anything that is not chemoinformatics goes here."""
# src/jazzy/helpers.py


def condense_atomic_map(atomic_map: dict) -> list:
    """Create condensed representation of the polar strength map.

    Polar strength map generated from `calculate_polar_strength_map()`.

    Args:
    atomic_map: polar strength map for all atoms in the system.

    """
    condensed_map = list()
    for idx, props in atomic_map.items():
        props["idx"] = idx
        condensed_map.append(props)
    return condensed_map


def convert_map_to_tuples(atomic_map: dict) -> list:
    """Create tuple representation of polar strength map.

    Tuple representation of the polar strength map generated from
    `calculate_polar_strength_map()`. Simple example:

    Args:
    atomic_map: polar strength map for all atoms in the system.

    Returns:
    List of tuples, where elements are atom indices and tuples are tuples of
    properties.

    """
    atomic_map_values = atomic_map.values()
    if len(atomic_map_values) == 0:
        raise IndexError("The atomic map must have length greater than zero.")

    tuple_map = list()
    for props in atomic_map_values:
        tuple_map.append(tuple((k, v) for k, v in props.items()))
    return tuple_map


def sum_atomic_map(atomic_map: dict) -> dict:
    """Optimized for summing the values within a polar strength map."""
    atomic_map_values = atomic_map.values()
    if len(atomic_map_values) == 0:
        raise IndexError("The atomic map must have length greater than zero.")

    # sum and round the countable values
    dct = dict()
    countable = ["sdc", "sdx", "sa"]
    for k in countable:
        dct[k] = sum(d[k] for d in atomic_map_values if d)
    return {key: round(dct[key], 4) for key in dct}

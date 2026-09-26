def gc_content(seq: str) -> float:
    """GC content as a percentage (0-100), rounded to 2 decimals.

    Case-insensitive; empty string and N-only sequences return 0.0.
    """
    if not isinstance(seq, str):
        raise TypeError("seq must be a single string")
    if len(seq) == 0:
        return 0.0
    s = seq.upper()
    return round(sum(base in "GC" for base in s) / len(s) * 100, 2)

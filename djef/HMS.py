# HMS module: converts string time input to float

def str_HMS(s):
    """Convert a string "H:M:S" to a tuple of numbers (H, M, S).
    H and M must be integers, S may be a float.
    """
    L = s.split(":")
    if len(L) < 3:
        raise ValueError("Too few fields in H:M:S string.")
    elif len(L) > 3:
        raise ValueError("Too many fields in H:M:S string.")
    H = int(L[0])
    M = int(L[1])
    S = float(L[2])
    return (H, M, S)

def HMS_Seconds(hours, minutes, seconds):
    """Convert hours minutes seconds to seconds."""
    return hours*60*60 + minutes*60 + seconds

def HMS_Hours(hours, minutes, seconds):
    """Convert hours minutes seconds to hours."""
    return hours + minutes/60.0 + seconds/(60.0**2)

def Hours_HMS(h):
    """Convert time t in hours to hours minutes seconds."""
    hours = int(t)
    t = (t - hours)*60
    minutes = int(t)
    seconds = (t - minutes)*60
    return (hours, minutes, seconds)

def Seconds_HMS(h):
    """Convert time t in seconds to hours minutes seconds."""
    hours, t = divmod(t, 60*60)
    minutes, seconds = divmod(t, 60)
    return (hours, minutes, seconds)

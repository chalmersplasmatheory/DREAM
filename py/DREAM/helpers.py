# Various helpful routines


def merge_dicts(old, new):
    """
    Update entries in the dictionary 'old' with the corresponding
    values in the dictionary 'new'. This function operates recursively.
    """
    d = old.copy()
    for k in new:
        if type(new[k]) == dict:
            if k not in d:
                d[k] = {}

            d[k] = merge_dicts(d[k], new[k])
        else:
            d[k] = new[k]

    return d


def safeTeXstring(s):
    mappings = {
        '#': r'\#',
        '_': r'\_'
    }

    for key, val in mappings.items():
        s = s.replace(key, val)
    
    return s

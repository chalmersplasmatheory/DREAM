# Various helpful routines


def safeTeXstring(s):
    mappings = {
        '#': r'\#',
        '_': r'\_'
    }

    for key, val in mappings.items():
        s = s.replace(key, val)
    
    return s

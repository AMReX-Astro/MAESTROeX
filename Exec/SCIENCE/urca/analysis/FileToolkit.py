from itertools import zip_longest

def compare_mask(s1, s2):
    """Compare string s1 to string s2, return elementwise equality as boolean mask."""

    return [a==b for a, b in zip_longest(str(s1), str(s2))]

def first_diff_substring(s, sref):
    """Return the first substring of string s which differs from the reference string sref.
    Return the empty string if the two strings are identical."""

    mask = compare_mask(s, sref)
    set_ifirst = False
    ifirst = 0
    ilast  = len(s)-1
    for i, x in enumerate(mask):
        if not x and not set_ifirst:
            ifirst = i
            set_ifirst = True
        elif x and set_ifirst:
            ilast = i-1
            break

    if set_ifirst and ifirst <= len(s)-1:
        return s[ifirst:ilast+1]
    else:
        return ''

def longest_diff_substring_list(s, sref_list):
    """Return the longest substring of string s which differs from all strings in sref_list.
    Return the empty string if no diffs can be found."""

    maxlen = 0
    maxdff = ''
    for sref in sref_list:
        sdiff = first_diff_substring(s, sref)
        if sdiff:
            if len(sdiff) > maxlen:
                maxdff = sdiff
                maxlen = len(sdiff)
    return maxdff

def get_sort_key(s, sref_list):
    """Return an integer sorting key for string s given the list of strings sref_list."""

    sdiff = longest_diff_substring_list(s, sref_list)
    if sdiff:
        return int(sdiff)
    else:
        return 0

def sort_input_filenames(filenames):
    """If sorting is desired, sort the input file list.
    Determine the substring to interpret as an integer by comparing file names."""

    if len(filenames) > 1:
        filenames = sorted(filenames, key=lambda x: get_sort_key(x, filenames))
    return filenames

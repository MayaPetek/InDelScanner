"""
A collection of random string utility functions. Mostly trivial.
"""


def getChar(target, index):
    """
    Get a char from a string, or "-" if `index` is out of range.
    """

    if index < 0 or index >= len(target):
        return "-";

    return target[index]

def getChars(target, index, n):
    """
    Get a sequence of `n` chars from `target` at `index`, safely, via `getChar`
    """

    out = ""
    for i in range(0, n):
        out += getChar(target, index + i);

    return out;

def setChar(targetStr, char, i):
    """
    Set the character at index `i` in `targetStr` to `char`, returning the new string.
    """

    return targetStr[:i] + char + targetStr[i + 1:]

def stealLetter(fromStr, toStr, i):
    """
    Steal a letter from `fromStr` and write it to `toStr` at `i`.
    """

    return setChar(toStr, fromStr[i], i);

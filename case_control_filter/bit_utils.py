'''Convenience methods for bit setting/checking.'''


def first_n_bits_set(x, n):
    '''
    Return True if x is equal to value of an int with first n bits set.

    e.g.
        all_bits_set(1, 1)  # evaluates to True
        all_bits_set(2, 2)  # evaluates to False
        all_bits_set(3, 2)  # evaluates to True
        all_bits_set(3, 1)  # evaluates to False
        all_bits_set(7, 3)  # evaluates to True
    '''
    return x == (1 << n) - 1


def set_first_bits(n):
    ''' Set first n bits to 1 '''
    return (1 << n) - 1


def set_bits_in_range(x, y):
    '''
    Set bits starting at 0-based index x up until y

    e.g.
        set_bits_in_range(0, 2)  # returns 3 (0b11)
        set_bits_in_range(1, 3)  # returns 6 (0b110)
        set_bits_in_range(2, 6)  # returns 60 (0b111100)

    '''
    return ((1 << (x)) - 1) ^ ((1 << (y)) - 1)


def flip_bits(x, n):
    '''Flip the first n bits of x'''
    return x ^ set_first_bits(n)


def highest_set_bit(x):
    ''' Return 0-based index of highest set bit in x'''
    if x == 0:
        return 0
    hsb = 0
    x = int(x / 2)
    while x > 0:
        x = int(x / 2)
        hsb += 1
    return hsb


def flag_consensus(flags):
    ''' For list of flags return bitwise & for all. '''
    f = flags[0]
    for x in flags[1:]:
        f &= x
    return f

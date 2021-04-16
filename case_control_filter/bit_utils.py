def all_bits_set(x, nbits):
    ''' Return True if first nbits bits are set. '''
    return all(x & 1 << y for y in range(nbits))


def set_bits_in_range(nbits):
    ''' Set first n bits to 1 '''
    x = 0
    for i in range(nbits):
        x |= (1 << i)
    return x


def flip_bits(x, nbits):
    '''Flip the first nbits bits of x'''
    return x ^ set_bits_in_range(nbits)


def highest_set_bit(x):
    ''' Return index of highest set bit in x'''
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

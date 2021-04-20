from nose.tools import *
from case_control_filter.bit_utils import first_n_bits_set, flag_consensus
from case_control_filter.bit_utils import flip_bits, highest_set_bit
from case_control_filter.bit_utils import set_bits_in_range, set_first_bits


def test_first_n_bits_set():
    assert(first_n_bits_set(0b1, 1))
    assert_false(first_n_bits_set(0b0, 1))
    assert(first_n_bits_set(0b1111, 4))
    assert_false(first_n_bits_set(0b111, 4))
    assert_false(first_n_bits_set(0b1011, 4))


def test_flag_consensus():
    assert_equal(flag_consensus([0b0, 0b11, 0b111]), 0)
    assert_equal(flag_consensus([0b1, 0b11, 0b111]), 1)
    assert_equal(flag_consensus([0b101, 0b1110, 0b111]), 0b100)


def test_flip_bits():
    assert_equal(flip_bits(0, 4), 0b1111)
    assert_equal(flip_bits(0b1010, 4), 0b101)
    assert_equal(flip_bits(0b101, 4), 0b1010)


def test_highest_set_bit():
    assert_equal(highest_set_bit(0), 0)
    assert_equal(highest_set_bit(1), 0)
    assert_equal(highest_set_bit(2), 1)
    assert_equal(highest_set_bit(0b100), 2)


def test_set_bits_in_range():
    assert_equal(set_bits_in_range(0, 3), 0b111)
    assert_equal(set_bits_in_range(2, 4), 0b1100)


def test_set_first_bits():
    assert_equal(set_first_bits(0), 0)
    assert_equal(set_first_bits(1), 1)
    assert_equal(set_first_bits(2), 3)
    assert_equal(set_first_bits(3), 0b111)


if __name__ == '__main__':
    import nose
    nose.run(defaultTest=__name__)

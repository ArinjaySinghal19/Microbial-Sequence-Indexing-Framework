#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 13:02:29 2023
@author: sai
"""

import mmh3
from bitarray import bitarray
import helper


class BloomFilter:
    """
    Bloom Filter class with add, find, intersection and union operations. Uses murmurhash
    functions for the operations.
    
    """
    def __init__(self, size: int, num_hashes: int, seed=42) -> None:
        self.size = size  # bit array size
        self.num_hashes = num_hashes  # number of hash functions
        self.seed = seed  # seed for mmh3
        self.bit_array = bitarray(self.size)  # bit array
        self.bit_array.setall(0)  # initialize all bits to 0
        self.salt=113

    def add(self, item: str) -> None:
        """
        Receives an item, string representation of 64-bit hash value
        and add its to the Bloom filter
        Parameters
        ----------
        item : str
            64-bit hash value.

        Returns
        -------
        None

        """
        for i in range(self.num_hashes):
            # create digests for given item with different seeds to each of mmh3 functions
            #digest = mmh3.hash64(item, self.seed+i, signed=False) % self.size
            digest = helper.get64bit(
                mmh3.hash64(item, self.seed+(self.salt*i), signed=False)) % self.size
            # set the bits in bit_array
            self.bit_array[digest] = True

    def find(self, item: str) -> bool:
        """
        Receives an item, string representation of 64-bit hash value
        and returns True, if the item is present in the Bloom filter
        else returns False.

        Parameters
        ----------
        item : str
            64-bit hash value.

        Returns
        -------
        bool
            True, if item is present, else False.
        """
        for i in range(self.num_hashes):
            #digest = mmh3.hash64(item, self.seed+i, signed=False) % self.size
            digest = helper.get64bit(
                mmh3.hash64(item, self.seed+(self.salt*i), signed=False)) % self.size
            if self.bit_array[digest] == False:
                return False
        return True

    def intersection(self, other):
        """
        Performs intersection of the Bloom Filter with the other Bloom Filter and
        returns a new Bloom Filter which represents the intersection of the two.

        Parameters
        ----------
        other : BloomFilter

        Returns
        -------
        res : BloomFilter (bit-wise intersection of two bloom filters)

        """
        res = BloomFilter(self.size, self.num_hashes)
        res.bit_array = self.bit_array & other.bit_array
        return res

    def union(self, other):
        """
        Performs union of the Bloom Filter with the other Bloom Filter and
        returns a new Bloom Filter which represents the union of the two.

        Parameters
        ----------
        other : BloomFilter

        Returns
        -------
        res : BloomFilter (bit-wise union)

        """
        res = BloomFilter(self.size, self.num_hashes)
        res.bit_array = self.bit_array | other.bit_array
        return res

    def add_all(self, items: set) -> None:
        """
        Receives a set of items (64-bit hash values) and
        add these items to Bloom Filter

        Parameters
        ----------
        items : set

        Returns
        -------
        None

        """
        for j in items:
            self.add(str(j))

    def count(self) -> (int, int):
        """
        Returns the count of bits that are turned on (1) and turned off (0) as a tuple
        in the Bloom Filter

        Returns
        -------
        (int, int) - (bit value 1, bit value 0)

        """
        ones = self.bit_array.count(1)
        return ones, self.size-ones

    def __or__(self, other) -> None:
        self.bit_array|=other.bit_array

    def __contains__(self, item: int) -> bool:
        return self.find(item)

    def __intersection__(self, other) -> None:
        self.bit_array = self.bit_array & other.bit_array

    def __union__(self, other) -> None:
        self.bit_array = self.bit_array | other.bit_array

    def reset(self) -> None:
        """
        Clears/Resets all bits in the Bloom Filter

        Returns
        -------
        None

        """
        self.bit_array.setall(0)  # initialize all bits to 0

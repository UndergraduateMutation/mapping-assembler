# Algorithms for mapping assembly

This project contains 3 algorithms for mapping reads to a reference and assembling them.

1. Naive approach
    - Take a seed from a read
    - Find its position in the reference
    - Expand from both sides until the mismatch count exceeds a threshold `t`
    - Save the index of the read, position in reference, the start and end indices in the read with max `2t` mismatches
    - Repeat for all reads
    - Sort by position in reference
    - Assemble genome

2. Knuth Morris Pratt
    - Find all the occurences of a read in reference
    - Sort by position in reference
    - Assemble genome

3. Suffix Array
    - Build a suffix array using DC3
    - Find all the occurences of a read in reference by binary searching the suffix array
    - Sort by position in reference
    - Assemble genome

## Quick start

To compile: 

```console
make
```

To run: 

```console
./main <path/to/reads> <path/to/reference> <algorithm>
```

The algorithm argument can be of the value `sa` for suffix array, `naive` for the naive approach, and `kmp` for KMP.

By default the suffix array solution will be used.

## References

Juha Kärkkäinen, J., Sanders, P., & Burkhardt, S. (n.d.). Linear work suffix array construction - University of Helsinki. https://www.cs.helsinki.fi/u/tpkarkka/publications/jacm05-revised.pdf

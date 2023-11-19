#pragma once

#include <string>
#include <vector>

// https://www.cs.helsinki.fi/u/tpkarkka/publications/jacm05-revised.pdf
inline bool leq(int a1, int a2, int b1, int b2);
inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3);
void radixPass(int *a, int *b, int *r, int n, int K);
void suffixArray(int *T, int *SA, int n, int K);
std::vector<int> build_suffix_array(const std::string &s);
void find_occurences(const std::string& pattern, const std::string& text, const std::vector<int>& suffix_array, std::vector<int>& result);
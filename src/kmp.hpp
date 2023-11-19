#pragma once

#include <string>
#include <vector>

void compute_lps(const std::string& pattern, int n, std::vector<int>& lps);
void kmp(const std::string& pattern, const std::string& text, std::vector<int>& occurences);
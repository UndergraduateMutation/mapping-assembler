#include "kmp.hpp"

void compute_lps(const std::string& pattern, int n, std::vector<int>& lps) {
	int len = 0;
	lps[0] = 0;
	int i = 1;
	while (i < n) {
		if (pattern[i] == pattern[len]) {
			len++;
			lps[i] = len;
			i++;
		} else {
			if (len != 0) {
				len = lps[len - 1];
			} else {
				lps[i] = 0;
				i++;
			}
		}
	}
}

void kmp(const std::string& pattern, const std::string& text, std::vector<int>& occurences) {
	int m = pattern.length();
	int n = text.length();

	std::vector<int> lps(m, 0);

	compute_lps(pattern, m, lps);

	int i = 0;
	int j = 0;
	while ((n - i) >= (m - j)) {
		if (pattern[j] == text[i]) {
			j++;
			i++;
		}

		if (j == m) {
			occurences.push_back(i - j);
			j = lps[j - 1];
		} else if (i < n && pattern[j] != text[i]) {
			if (j != 0)
				j = lps[j - 1];
			else
				i = i + 1;
		}
	}
}
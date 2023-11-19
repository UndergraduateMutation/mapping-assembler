#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <tuple>
#include <unordered_map>

enum alg {
	naive, sa, knuth
};

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

void sort_characters(const std::string& text, std::vector<int>& order) {
	std::vector<int> count(5, 0);
	std::map<char, int> charToIdx;
	charToIdx['$'] = 0;
	charToIdx['a'] = 1;
	charToIdx['c'] = 2;
	charToIdx['g'] = 3;
	charToIdx['t'] = 4;
	for (size_t i = 0; i < text.size(); i++) {
		count[charToIdx[text[i]]]++;
	}
	for (size_t i = 1; i < count.size(); i++) {
		count[i] += count[i - 1];
	}
	for (int i = text.size() - 1; i >= 0; i--) {
		int idx = charToIdx[text[i]];
		count[idx] -= 1;
		order[count[idx]] = i;
	}
}

void compute_labels(const std::string& text, const std::vector<int>& order, std::vector<int>& label) {
	label[order[0]] = 0;
	for (size_t i = 1; i < text.length(); i++) {
		if (text[order[i]] != text[order[i - 1]]) {
			label[order[i]] = label[order[i - 1]] + 1;
		} else {
			label[order[i]] = label[order[i - 1]];
		}
	}
}

void sorted_doubled(const std::string& text, int L, std::vector<int>& order, const std::vector<int>& label) {
	std::vector<int> newOrder(text.size());
	std::vector<int> count(text.size(), 0);
	for (size_t i = 0; i < text.length(); i++) {
		count[label[i]]++;
	}
	for (size_t i = 1; i < count.size(); i++) {
		count[i] += count[i - 1];
	}
	for (int i = text.size() - 1; i >= 0; i--) {
		int start = (order[i] - L + text.size()) % text.size();
		int cl = label[start];
		count[cl] -= 1;
		newOrder[count[cl]] = start;
	}
	order = newOrder;
}

void update_labels(const std::vector<int>& newOrder, std::vector<int>& label, int L) {
	std::vector<int> newlabel(newOrder.size());
	newlabel[newOrder[0]] = 0;
	for (size_t i = 1; i < newOrder.size(); i++) {
		int cur = newOrder[i];
		int prev = newOrder[i - 1];
		int mid = cur + L;
		int midPrev = (prev + L) % newOrder.size();
		if (label[cur] != label[prev] || label[mid] != label[midPrev]) {
			newlabel[cur] = newlabel[prev] + 1;
		} else {
			newlabel[cur] = newlabel[prev];
		}
	}
	label = newlabel;
}

std::vector<int> build_sa(const std::string& text) {
	std::vector<int> result;
	std::vector<int> order(text.size());
	std::vector<int> label(text.size());
	sort_characters(text, order);
	compute_labels(text, order, label);
	size_t L = 1;
	while (L < text.length()) {
		sorted_doubled(text, L, order, label);
		update_labels(order, label, L);
		L *= 2;
	}
	result = order;
	return result;
}

void find_occurences(const std::string& pattern, const std::string& text, const std::vector<int>& suffix_array, std::vector<int>& result) {
	int minIndex = 0;
	int maxIndex = text.size();
	while(minIndex < maxIndex) {
		int midIndex = (minIndex + maxIndex) / 2;
		int flag = text.compare(suffix_array[midIndex], std::min(text.size() - suffix_array[midIndex], pattern.size()), pattern);
		if (flag < 0) {
			minIndex = midIndex + 1;
		} else {
			maxIndex = midIndex;
		}
	}

	int start = minIndex;
	maxIndex = text.size();
	while(minIndex < maxIndex) {
		int midIndex = (minIndex + maxIndex) / 2;
		int flag = text.compare(suffix_array[midIndex], std::min(text.size() - suffix_array[midIndex], pattern.size()), pattern);
		if (flag > 0) {
			maxIndex = midIndex;
		} else {
			minIndex = midIndex + 1;
		}
	}

	int end = maxIndex;
	if (start <= end) {
		for (int i = start; i < end; i++) {
			result.push_back(suffix_array[i]);
		}
	}
}

std::string get_consensus(const std::vector<std::string>& reads) {
    if (reads.empty()) {
        return "";
    }

    int read_length = reads[0].size();
    std::string consensus(read_length, ' ');

    for (int i = 0; i < read_length; i++) {
        std::unordered_map<char, int> counts;

        for (const std::string& read : reads) {
            char base = read[i];
            if (counts.find(base) == counts.end()) {
                counts[base] = 1;
            } else {
                counts[base]++;
            }
        }

        char most_freq_base = 'N';
        int max_count = 0;

        for (const auto& pair : counts) {
            if (pair.second > max_count) {
                max_count = pair.second;
                most_freq_base = pair.first;
            }
        }

        consensus[i] = most_freq_base;
    }

    return consensus;
}

void assemble_sa(std::string& reference, std::vector<std::string>& reads, std::ofstream& out){
	std::vector<int> suffix_array = build_sa(reference);
	int index = 0;
	std::vector<std::pair<int, int>> assembledReads;
    for (const auto& read : reads) {
		std::vector<int> occurrences;
		find_occurences(read, reference, suffix_array, occurrences);
		for (size_t j = 0; j < occurrences.size(); j++) {
			assembledReads.push_back({occurrences[j], index});
		}
		index++;
    }

	std::sort(assembledReads.begin(), assembledReads.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

	std::stringstream assembly;
	int cursor = assembledReads[0].first;
	for(const auto& p: assembledReads){
		if(cursor > p.first){
			size_t overlap = cursor - p.first;
			std::string as = assembly.str();
			std::string consensus = get_consensus({as.substr(as.length() - overlap, as.length()), reads[p.second].substr(0, overlap)});
			as = as.substr(0, as.length() - overlap);
			assembly.str("");
			assembly.clear();
			assembly << as;
			assembly << consensus;
			assembly << reads[p.second].substr(overlap, reads[p.second].length());
			cursor += reads[p.second].length() - overlap;
		} else {
			assembly << reads[p.second];
			cursor += reads[p.second].length();
		}
	}

	out << assembly.str();
}

void assemble_kmp(std::string& reference, std::vector<std::string>& reads, std::ofstream& out) {
	int index = 0;
	std::vector<std::pair<int, int>> assembledReads;
    for (const auto& read : reads) {
        std::vector<int> occurrences;
		kmp(read, reference, occurrences);
		for (size_t j = 0; j < occurrences.size(); j++) {
			assembledReads.push_back({occurrences[j], index});
		}
		index++;
    }

	std::sort(assembledReads.begin(), assembledReads.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

	std::stringstream assembly;
	int cursor = assembledReads[0].first;
	for(const auto& p: assembledReads){
		if(cursor > p.first){
			size_t overlap = cursor - p.first;
			std::string as = assembly.str();
			std::string consensus = get_consensus({as.substr(as.length() - overlap, as.length()), reads[p.second].substr(0, overlap)});
			as = as.substr(0, as.length() - overlap);
			assembly.str("");
			assembly.clear();
			assembly << as;
			assembly << consensus;
			assembly << reads[p.second].substr(overlap, reads[p.second].length());
			cursor += reads[p.second].length() - overlap;
		} else {
			assembly << reads[p.second];
			cursor += reads[p.second].length();
		}
	}

	out << assembly.str();
}

void assemble_naive(std::string& reference, std::vector<std::string>& reads, std::ofstream& out) {
	int index = 0;
	size_t max_mismatches = 5;
	size_t seed_len = 5;
	std::vector<std::tuple<int, int, int, int>> matching_reads;
	for (const auto& read : reads) {
        int rand_index = rand() % (seed_len + 1);
        size_t i = reference.find(read.substr(rand_index, seed_len));
        if (i != std::string::npos) {
            size_t count = 0, read_end, ref_end, read_start, ref_start, mismatch;
            for(read_start = rand_index - 1, ref_start = i - 1; 
                read_start > 0 && ref_start > 0; 
                read_start--, ref_start--) {
                if(count > max_mismatches) break;

                if(read[read_start]!=reference[ref_start]) count++;
            }
            mismatch = count;
            count = 0;
            for(read_end = rand_index + seed_len, ref_end = i + seed_len; 
                read_end < read.length() && ref_end < reference.length(); 
                read_end++, ref_end++) {
                if(count > max_mismatches) break;

                if(read[read_end]!=reference[ref_end]) count++;
            }
            mismatch+=count;

			matching_reads.push_back({ref_start, read_start, read_end, index});
			index++;
        }
    }

	std::sort(matching_reads.begin(), matching_reads.end(), [](const auto& a, const auto& b) {
        return std::get<0>(a) < std::get<0>(b);
    });

	std::stringstream assembly;
	int cursor = std::get<0>(matching_reads[0]);
	for(const auto& p: matching_reads){
		std::string matching_string = reads[std::get<3>(p)].substr(std::get<1>(p), std::get<2>(p));
		if(cursor > std::get<0>(p)){
			size_t overlap = cursor - std::get<0>(p);
			std::string as = assembly.str();
			std::string consensus = get_consensus({as.substr(as.length() - overlap, as.length()), matching_string.substr(0, overlap)});
			as = as.substr(0, as.length() - overlap);
			assembly.str("");
			assembly.clear();
			assembly << as;
			assembly << consensus;
			assembly << matching_string.substr(overlap, matching_string.length());
			cursor += matching_string.length() - overlap;
		} else {
			assembly << matching_string;
			cursor += matching_string.length();
		}
	}

	out << assembly.str();
}

void assemble(std::string& reference, std::vector<std::string>& reads, std::ofstream& out, const alg algo) {
	switch(algo) {
		case alg::sa:
			assemble_sa(reference, reads, out);
			return;
		case alg::knuth:
			assemble_kmp(reference, reads, out);
			return;
		case alg::naive:
			assemble_naive(reference, reads, out);
			return;
		default: 
			return;
	}
}

int main() {
    std::ifstream in_ref("reference.txt");
    std::ifstream in_reads("reads.txt");
    std::ofstream out("output.txt");
    std::stringstream ref;

    std::vector<std::string> reads;

    std::string line;
    while (in_reads >> line) {
        reads.push_back(line);
    }

    while (std::getline(in_ref, line)) {
        ref << line;
    }

    std::string reference = ref.str();
    ref.str(""); 
    ref.clear();

	assemble(reference, reads, out, alg::naive);
   
    in_ref.close();
    in_reads.close();

    return 0;
}

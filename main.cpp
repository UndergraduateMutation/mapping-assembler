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

// https://www.cs.helsinki.fi/u/tpkarkka/publications/jacm05-revised.pdf
inline bool leq(int a1, int a2, int b1, int b2) {
    return (a1 < b1 || (a1 == b1 && a2 <= b2));
}

// lexicographic order for triples
inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3) {
    return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
}

// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
static void radixPass(int *a, int *b, int *r, int n, int K) {
    // count occurrences
    int *c = new int[K + 1];  // counter array
    for (int i = 0; i <= K; i++)
        c[i] = 0;  // reset counters
    for (int i = 0; i < n; i++)
        c[r[a[i]]]++;                      // count occurrences
    for (int i = 0, sum = 0; i <= K; i++)  // exclusive prefix sums
    {
        int t = c[i];
        c[i] = sum;
        sum += t;
    }
    for (int i = 0; i < n; i++)
        b[c[r[a[i]]]++] = a[i];  // sort
    delete[] c;
}

// find the suffix array SA of T[0..n-1] in {1..K}^n
// require T[n]=T[n+1]=T[n+2]=0, n>=2
void suffixArray(int *T, int *SA, int n, int K) {
    int n0 = (n + 2) / 3, n1 = (n + 1) / 3, n2 = n / 3, n02 = n0 + n2;
    int *R = new int[n02 + 3];
    R[n02] = R[n02 + 1] = R[n02 + 2] = 0;
    int *SA12 = new int[n02 + 3];
    SA12[n02] = SA12[n02 + 1] = SA12[n02 + 2] = 0;
    int *R0 = new int[n0];
    int *SA0 = new int[n0];

    //******* Step 0: Construct sample ********
    // generate positions of mod 1 and mod 2 suffixes
    // the "+(n0-n1)" adds a dummy mod 1 suffix if n%3 == 1
    for (int i = 0, j = 0; i < n + (n0 - n1); i++)
        if (i % 3 != 0)
            R[j++] = i;

    //******* Step 1: Sort sample suffixes ********
    // lsb radix sort the mod 1 and mod 2 triples
    radixPass(R, SA12, T + 2, n02, K);
    radixPass(SA12, R, T + 1, n02, K);
    radixPass(R, SA12, T, n02, K);

    // find lexicographic names of triples and
    // write them to correct places in R
    int name = 0, c0 = -1, c1 = -1, c2 = -1;
    for (int i = 0; i < n02; i++) {
        if (T[SA12[i]] != c0 || T[SA12[i] + 1] != c1 || T[SA12[i] + 2] != c2) {
            name++;
            c0 = T[SA12[i]];
            c1 = T[SA12[i] + 1];
            c2 = T[SA12[i] + 2];
        }
        if (SA12[i] % 3 == 1) {
            R[SA12[i] / 3] = name;
        }  // write to R1
        else {
            R[SA12[i] / 3 + n0] = name;
        }  // write to R2
    }

    // recurse if names are not yet unique
    if (name < n02) {
        suffixArray(R, SA12, n02, name);
        // store unique names in R using the suffix array
        for (int i = 0; i < n02; i++)
            R[SA12[i]] = i + 1;
    } else  // generate the suffix array of R directly
        for (int i = 0; i < n02; i++)
            SA12[R[i] - 1] = i;
    //******* Step 2: Sort nonsample suffixes ********
    // stably sort the mod 0 suffixes from SA12 by their first character
    for (int i = 0, j = 0; i < n02; i++)
        if (SA12[i] < n0)
            R0[j++] = 3 * SA12[i];
    radixPass(R0, SA0, T, n0, K);

    //******* Step 3: Merge ********
    // merge sorted SA0 suffixes and sorted SA12 suffixes
    for (int p = 0, t = n0 - n1, k = 0; k < n; k++) {
#define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
        int i = GetI();     // pos of current offset 12 suffix
        int j = SA0[p];     // pos of current offset 0 suffix
        if (SA12[t] < n0 ?  // different compares for mod 1 and mod 2 suffixes
                leq(T[i], R[SA12[t] + n0], T[j], R[j / 3])
                         : leq(T[i], T[i + 1], R[SA12[t] - n0 + 1], T[j], T[j + 1],
                               R[j / 3 + n0])) {  // suffix from SA12 is smaller
            SA[k] = i;
            t++;
            if (t == n02)  // done --- only SA0 suffixes left
                for (k++; p < n0; p++, k++)
                    SA[k] = SA0[p];
        } else {  // suffix from SA0 is smaller
            SA[k] = j;
            p++;
            if (p == n0)  // done --- only SA12 suffixes left
                for (k++; t < n02; t++, k++)
                    SA[k] = GetI();
        }
    }
    delete[] R;
    delete[] SA12;
    delete[] SA0;
    delete[] R0;
}

std::vector<int> build_suffix_array(const std::string &s) {
    int n = s.size();
    if (n == 0)
        return {};
    if (n == 1)
        return {0};
    std::vector<int> t(n + 3);
    for (int i = 0; i < n; i++)
        t[i] = s[i];
    std::vector<int> sa(n);
    suffixArray(t.data(), sa.data(), n, 256);
    return sa;
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
	std::vector<int> suffix_array = build_suffix_array(reference);
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

int main(int argc, char *argv[]) {
	if(argc != 4){
		std::cout << "USAGE: " << argv[0] << " <path/to/reads> <path/to/reference> <algorithm>\n";
		std::cout << "ALGORITHMS: sa, naive, kmp\n";
		return 1;
	}

	std::string reads_path(argv[1]);
	std::string reference_path(argv[2]);
	std::string algorithm(argv[3]);

	alg algo = alg::sa;

	if(algorithm == "kmp") {
		algo = alg::knuth;
	} else if (algorithm == "naive"){
		algo = alg::naive;
	}

    std::ifstream in_ref(reference_path);
    std::ifstream in_reads(reads_path);
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

	assemble(reference, reads, out, algo);
   
    in_ref.close();
    in_reads.close();

    return 0;
}

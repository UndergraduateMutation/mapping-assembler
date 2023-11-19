#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <tuple>
#include <unordered_map>

#include "src/dc3.hpp"
#include "src/kmp.hpp"

enum alg {
	naive, sa, knuth
};

void assemble_sa(std::string& reference, std::vector<std::string>& reads, std::ofstream& out){
	std::vector<int> suffix_array = build_suffix_array(reference);
	int index = 0;
	std::vector<std::pair<int, int>> assembled_reads;
    for (const auto& read : reads) {
		std::vector<int> occurrences;
		find_occurences(read, reference, suffix_array, occurrences);
		for (size_t j = 0; j < occurrences.size(); j++) {
			assembled_reads.push_back({occurrences[j], index});
		}
		index++;
    }

	std::sort(assembled_reads.begin(), assembled_reads.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

	int cursor = assembled_reads[0].first;
	for(const auto& p: assembled_reads){
		if(cursor > p.first){
			size_t overlap = cursor - p.first;
			out << reads[p.second].substr(overlap, reads[p.second].length());
			cursor += reads[p.second].length() - overlap;
		} else {
			out << reads[p.second];
			cursor += reads[p.second].length();
		}
	}
}

void assemble_kmp(std::string& reference, std::vector<std::string>& reads, std::ofstream& out) {
	int index = 0;
	std::vector<std::pair<int, int>> assembled_reads;
    for (const auto& read : reads) {
        std::vector<int> occurrences;
		kmp(read, reference, occurrences);
		for (size_t j = 0; j < occurrences.size(); j++) {
			assembled_reads.push_back({occurrences[j], index});
		}
		index++;
    }

	std::sort(assembled_reads.begin(), assembled_reads.end(), [](const auto& a, const auto& b) {
        return a.first < b.first;
    });

	int cursor = assembled_reads[0].first;
	for(const auto& p: assembled_reads){
		if(cursor > p.first){
			size_t overlap = cursor - p.first;
			out << reads[p.second].substr(overlap, reads[p.second].length());
			cursor += reads[p.second].length() - overlap;
		} else {
			out << reads[p.second];
			cursor += reads[p.second].length();
		}
	}
}

void assemble_naive(std::string& reference, std::vector<std::string>& reads, std::ofstream& out) {
	int index = 0;
	size_t max_mismatches = 5;
	size_t seed_len = 5;
	std::vector<std::tuple<int, int, int, int>> matching_reads;
	for (const auto& read : reads) {
        int rand_index = rand() % (read.length() - seed_len);
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

	int cursor = std::get<0>(matching_reads[0]);
	for(const auto& p: matching_reads){
		std::string matching_string = reads[std::get<3>(p)].substr(std::get<1>(p), std::get<2>(p));
		if(cursor > std::get<0>(p)){
			size_t overlap = cursor - std::get<0>(p);
			out << matching_string.substr(overlap, matching_string.length());
			cursor += matching_string.length() - overlap;
		} else {
			out << matching_string;
			cursor += matching_string.length();
		}
	}
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

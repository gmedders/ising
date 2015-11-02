//example use case with uniformity test
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <string>

#include "random_selector.h"

typedef std::unordered_map<int, unsigned> histogram_t;

void run_testcase(std::vector<int>& source_container,
		  const unsigned samples = 10000u) {
    histogram_t histogram;
    random_selector<> selector{};
    for (auto i = samples; i > 0; --i) {
	++(histogram[selector(source_container)]);
    }

    std::cout << "Selection histogram from " << samples << " samples:\n";
    for (auto i : histogram) {
	std::cout << i.first << ':' << ' ' << i.second << '\n';
    }
}

//just some silly preamble
int main(int argc, const char* argv[]) {
    std::vector<int> source_container;
    source_container.push_back(2);
    source_container.push_back(43);
    source_container.push_back(443);

    run_testcase(source_container);
    return 0;
}

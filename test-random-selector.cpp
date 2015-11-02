//example use case with uniformity test
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iterator>
#include <string>

#include "random_selector.h"

typedef std::vector<std::string> vec_t;
typedef std::unordered_map<typename vec_t::value_type, unsigned> histogram_t;

void run_testcase(vec_t& source_container, const unsigned samples = 10000u) {
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
    vec_t source_container;
    if (argc > 1) {
	source_container.assign(argv+1, argv+argc);
    } else {
	source_container = {"2", "43", "443"};
    }

    run_testcase(source_container);
    return 0;
}

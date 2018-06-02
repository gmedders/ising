//example use case with uniformity test
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

void run_testcase(std::vector<int>& source_container,
		  const unsigned samples = 10000u)
{
    int counter[3] = {0, 0, 0};

//    std::random_device rd;
//    std::mt19937 g(rd());
	 
    for (auto i = samples; i > 0; --i) {
	std::random_shuffle(source_container.begin(), source_container.end());
	//std::random_shuffle(source_container.begin(), source_container.end(), g);
	for(auto j = 0; j < source_container.size(); ++j)
#if 0
//	    std::cout << shuffled[j] << ' ';
//	std::cout << std::endl;
#else
	    if(source_container[j] == 443)
		++counter[j];
#endif
    }
    for(auto i = 0; i < 3; ++i)
	std::cout << counter[i] << ' ';
    std::cout << std::endl;
}

int main(int argc, const char* argv[]) {
    std::vector<int> source_container;
    source_container.push_back(2);
    source_container.push_back(43);
    source_container.push_back(443);

//  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//
//  auto engine = std::default_random_engine{};

    run_testcase(source_container);
    return 0;
}

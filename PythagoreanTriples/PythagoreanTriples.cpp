#include <iostream>
#include <random>
#include"PythagoreanTriples.hpp"


int main()
{
    std::cout << "Hello World!\n\r";
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int<int> ud(0, 100);
    std::cout << "random number: " << ud(mt);
    //tpp::generate_triples(2llu * 5llu * 13llu * 17llu * 65537llu*100003llu, mt);
    //tpp::generate_triples(4611686018427387902, mt);
    //tpp::factor(2llu*5llu*13llu*17llu* 65537llu * 100003llu, mt);
    auto factorization = tpp::factor(4611686890305753337, mt);
    for (const auto& it : factorization)
        std::cout << " " << std::get<0>(it) << "^" << std::get<1>(it);
}
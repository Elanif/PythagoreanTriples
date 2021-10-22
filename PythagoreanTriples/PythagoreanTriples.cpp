#include <iostream>
#include <random>
#include"PythagoreanTriples.hpp"


int main()
{
    std::cout << "Hello World!\n\r";
    std::random_device rd;
    std::mt19937 mt(rd());
    //tpp::generate_triples(2llu * 5llu * 13llu * 17llu * 65537llu*100003llu, mt);
    std::cout << "int 32 max" << std::numeric_limits<int_fast32_t>::max()<<"\n\r";
    std::cout << "int 64 max" << std::numeric_limits<int_fast64_t>::max()<<"\n\r";
    std::cout << tpp::ispmiller(2305843009213693951ull, mt) * 10 << "\n\r";
    auto factorization = tpp::factor(4611686890305753337ull, mt);// 4611686890305753337, mt);
    for (const auto& it : factorization)
        std::cout << " " << std::get<0>(it) << "^" << std::get<1>(it);

    //big prime 2305843009213693951
    //big product of 2 primes 4611686890305753337
}
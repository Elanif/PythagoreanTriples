enum class multiplication_type {
    default_multiplication,
    swap_multiplication,
    increasing_multiplication,
    decreasing_multiplication,
    positive_multiplication,
    negative_multiplication
} multiplication;
#include <iostream>
#include <random>
#include"PythagoreanTriples.hpp"

int main()
{
    multiplication = multiplication_type::default_multiplication;
    std::cout << "Hello World!\n\r";
    std::random_device rd;
    std::mt19937 mt(rd());
    //tpp::generate_triples(2llu * 5llu * 13llu * 17llu * 65537llu*100003llu, mt);
    std::cout << "int 32 max " << std::numeric_limits<int_fast32_t>::max()<<"\n\r";
    std::cout << "int 64 max " << std::numeric_limits<int_fast64_t>::max()<<"\n\r";
    auto factorization = tpp::factor(uint_fast64_t(4611686890305753337ull), mt);
    for (const auto& it : factorization)
        std::cout << std::get<0>(it) << "^" << std::get<1>(it)<<" ";
    std::cout << "\n\r";
    tpp::pp<int_fast64_t> five(2, 1);
    tpp::pp<int_fast64_t> thirteen(3, 2);
    tpp::pp<int_fast64_t> seventeen(1, 4);
    tpp::pp<int_fast64_t> twentynine(5, 2);
    tpp::pp<int_fast64_t> unity(1, 0);
    /*std::cout << unity << "*" << five << "="<<unity * five<<"; "<<five<<"*"<<unity<<"="<<five*unity;
    std::cout << "\n\r";
    std::cout << unity << "*" << thirteen << "=" << unity * thirteen << "; " << thirteen << "*" << unity << "=" << thirteen * unity;
    std::cout << "\n\r";
    std::cout << five << "*" << thirteen << "=" << five * thirteen << "; " << thirteen << "*" << five << "=" << thirteen * five;
    std::cout << "\n\r";*/
    std::cout << five << "*" << five << "=" << five * five << "; " << five << "*" << five << "*" << five << "=" << five * five * five;
    std::cout << "\n\r";
    tpp::pp<int_fast64_t> parity_1 = five * thirteen;
    tpp::pp<int_fast64_t> parity_2 = parity_1* parity_1;
    std::cout << parity_1 << "," << parity_2<<"," << parity_1*parity_2<< "\n\r";

    parity_1 = five * thirteen*seventeen;
    /*std::cout << parity_1 << "," << parity_1*seventeen << parity_1 * seventeen * twentynine << parity_1 * seventeen * five << "\n\r";*/
    std::cout << parity_1 << seventeen << five<< "\n\r";;
    std::cout<< parity_1 * seventeen * five << "\n\r";
    std::cout << std::gcd((parity_1 * seventeen * five).a, (parity_1 * seventeen * five).b) << "\n\r";

    if (five.compatible(five)) std::cout << "five compatible with five\n\r";
    std::cout << five % five * five<< "\n\r";
    tpp::pp<int_fast64_t> test1 = tpp::make_pp<int_fast64_t>(2, 1);
    tpp::pp<int_fast64_t> test2 = tpp::make_pp<int_fast64_t>(24, 7);
    std::cout << test1 * test2;

    //big prime 2305843009213693951
    //big product of 2 primes 4611686890305753337

    //(1,4)*(1,8) and (1,4)*(8,1) both give good results
}
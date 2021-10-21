#pragma once
#include"random.hpp"
#include<map>
#include<numeric>
#include<array>
#include<random>



namespace tpp {
	using Random = effolkronium::random_static;

	typedef std::tuple<int_fast64_t, int_fast64_t> pair;
	typedef std::tuple<int_fast64_t, int_fast64_t, int_fast64_t> triple;

	int_fast64_t ipow(int_fast64_t base, uint_fast8_t exp) { //https://gist.github.com/orlp/3551590
		static const uint_fast8_t highest_bit_set[] = {
			0, 1, 2, 2, 3, 3, 3, 3,
			4, 4, 4, 4, 4, 4, 4, 4,
			5, 5, 5, 5, 5, 5, 5, 5,
			5, 5, 5, 5, 5, 5, 5, 5,
			6, 6, 6, 6, 6, 6, 6, 6,
			6, 6, 6, 6, 6, 6, 6, 6,
			6, 6, 6, 6, 6, 6, 6, 6,
			6, 6, 6, 6, 6, 6, 6, 255, // anything past 63 is a guaranteed overflow with base > 1
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
			255, 255, 255, 255, 255, 255, 255, 255,
		};

		int64_t result = 1;

		switch (highest_bit_set[exp]) {
		case 255: // we use 255 as an overflow marker and return 0 on overflow/underflow
			if (base == 1) {
				return 1;
			}

			if (base == -1) {
				return 1 - 2 * (exp & 1);
			}

			return 0;
		case 6:
			if (exp & 1) result *= base;
			exp >>= 1;
			base *= base;
		case 5:
			if (exp & 1) result *= base;
			exp >>= 1;
			base *= base;
		case 4:
			if (exp & 1) result *= base;
			exp >>= 1;
			base *= base;
		case 3:
			if (exp & 1) result *= base;
			exp >>= 1;
			base *= base;
		case 2:
			if (exp & 1) result *= base;
			exp >>= 1;
			base *= base;
		case 1:
			if (exp & 1) result *= base;
		default:
			return result;
		}
	}

	int_fast64_t mult(int_fast64_t a, int_fast64_t b, int_fast64_t mod) {
		int_fast64_t result = 0;
		while (b) {
			if (b & 1)
				result = (result + a) % mod;
			a = (a + a) % mod;
			b >>= 1;
		}
		return result;
	}

	int_fast64_t f(int_fast64_t x, int_fast64_t c, int_fast64_t mod) {
		return (mult(x, x, mod) + c) % mod;
	}

	//https://cp-algorithms.com/algebra/factorization.html
	int_fast64_t brent(int_fast64_t n, int_fast64_t x0 = 2, int_fast64_t c = 1) {
		int_fast64_t x = x0;
		int_fast64_t g = 1;
		int_fast64_t q = 1;
		int_fast64_t xs, y;

		int_fast64_t m = 64; // 128;
		int_fast64_t l = 1;
		while (g == 1) {
			y = x;
			for (int_fast64_t i = 1; i < l; i++)
				x = f(x, c, n);
			int_fast64_t k = 0;
			while (k < l && g == 1) {
				xs = x;
				for (int_fast64_t i = 0; i < m && i < l - k; i++) {
					x = f(x, c, n);
					q = mult(q, abs(y - x), n);
				}
				g = std::gcd(q, n);
				k += m;
			}
			l *= 2;
		}
		if (g == n) {
			do {
				xs = f(xs, c, n);
				g = std::gcd(abs(xs - y), n);
			} while (g == 1);
		}
		return g;
	}

	//https://cp-algorithms.com/algebra/factorization.html
	int_fast64_t trial_division3(int_fast64_t n, std::map< int_fast64_t, int_fast64_t>& previous, int_fast64_t max) {
		auto add = [&](int_fast64_t d) {
			auto it = previous.find(d);
			if (it != previous.end())
				++it->second;
			else
				previous[d] = 1;
		};
		for (int_fast64_t d : {2, 3, 5}) {
			while (n % d == 0) {
				add(d);
				n /= d;
			}
		}
		static std::array<int_fast64_t, 8> increments = { 4, 2, 4, 2, 4, 6, 2, 6 };
		int_fast64_t i = 0;
		for (int_fast64_t d = 7; d * d <= std::min(n,max); d += increments[i++]) {
			while (n % d == 0) {
				add(d);
				n /= d;
			}
			if (i == 8)
				i = 0;
		}
		return n;
	}
	
	template<class mersenne_twister_engine>
	void factor_recursive(int_fast64_t n, std::map< int_fast64_t, int_fast64_t>& previous, mersenne_twister_engine& mt, int_fast64_t depth) {
		if (n == 1) return;
		std::uniform_int<int_fast64_t> ud(2, n - 1);
		auto x0 = ud(mt);
		auto c = ud(mt);
		int_fast64_t result = brent(n, x0, c);
		auto it = previous.find(result);
		if (it != previous.end()) //if it's in then either it's fucked or it's ok to do this
			it->second++;
		else {
			if (result == n || result == 1) {
				if (depth <= 0) {// give up
					previous[n] = 1;
				}
				else { //try again
					factor_recursive(n, previous, mt, depth - 1);
				}
			}
			else {
				factor_recursive(n / result, previous, mt, depth);
				factor_recursive(result, previous, mt, depth);
			}
		}
	}

	template<class mersenne_twister_engine>
	std::vector<pair> factor(int_fast64_t n, mersenne_twister_engine& mt) {
		std::vector<pair> result;
		std::map< int_fast64_t, int_fast64_t> factorization_map;
		//factorization_map[1] = 1;
		n = trial_division3(n, factorization_map, 100000);
		if (n>1)
			factor_recursive(n, factorization_map, mt, 1);
		result.assign(factorization_map.begin(), factorization_map.end());
		return result;
	};

	std::tuple<pair> euclidian_algorithm(int_fast64_t big, int_fast64_t small, int_fast64_t cutoff) {
		return pair(1, 0);
	};

	std::vector<triple> generate_triples(const std::vector<pair>& factorization) {
		int_fast64_t multiplier = 1;
		std::vector<pair> good_primes;
		for (auto const& it : factorization) {
			if (std::get<0>(it) % 4 == 1) {
				good_primes.push_back(it);
			}
			else {
				multiplier *= ipow(std::get<1>(it), std::get<0>(it));
			}
		}
		std::vector<triple> result;
		return result;
	}

	template<class mersenne_twister_engine>
	std::vector<triple> generate_triples(int_fast64_t n, mersenne_twister_engine& mt) {
		return generate_triples(factor(n, mt));
	}
}
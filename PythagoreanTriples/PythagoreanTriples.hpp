#pragma once
#include"random.hpp"
#include<map>
#include<numeric>
#include<array>
#include<random>

namespace tpp {

	using Random = effolkronium::random_static;

	template<typename int_type>
	using pair = std::tuple<int_type, int_type>;
	template<typename int_type>
	using triple = std::tuple<int_type, int_type, int_type>;

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

	template<typename int_type>
	int_type mul_mod(int_type a, int_type b, int_type mod) {
		int_type result = 0;
		while (b) {
			if (b & 1)
				result = (result + a) % mod;
			a = (a + a) % mod;
			b >>= 1;
		}
		return result;
	}

	template<typename int_type>
	int_type power_mod(int_type a, int_type b, int_type mod) {
		int_type result = 1;
		while (b) {
			if (b & 1) result = mul_mod(result, a, mod);
			if (b >>= 1) a = mul_mod(a, a, mod);
		}
		return result % mod;
	}

	template<typename int_type>
	int_type f(int_type x, int_type c, int_type mod) {
		return (mul_mod(x, x, mod) + c) % mod;
	}

	template<typename int_type>
	int_type abs(int_type x, int_type y) {
		if (x < y) return y - x;
		return x - y;
	}

	//https://cp-algorithms.com/algebra/factorization.html
	template<typename int_type>
	int_type trial_division3(int_type n, std::map<int_type, int_type>& previous, int_type max) {
		auto add = [&](int_type d) {
			auto it = previous.find(d);
			if (it != previous.end())
				++it->second;
			else
				previous[d] = 1;
		};
		for (int_type d : {2, 3, 5}) {
			while (n % d == 0) {
				add(d);
				n /= d;
			}
		}
		static std::array<int_type, 8> increments = { 4, 2, 4, 2, 4, 6, 2, 6 };
		int_type i = 0;
		for (int_type d = 7; d * d <= std::min(n, max); d += increments[i++]) {
			while (n % d == 0) {
				add(d);
				n /= d;
			}
			if (i == 8)
				i = 0;
		}
		return n;
	}

	template<typename int_type, class mersenne_twister_engine>
	bool ispmiller(int_type p, mersenne_twister_engine& mt) { // O(30*logp)
		if (p < 2) return false;
		if (p == 2) return true;
		if (p % 2 == 0) return false;
		std::uniform_int_distribution<int_type> ud(3, p - 1);
		int_type s = p - 1;
		unsigned long trailing_zeroes = 0;
		_BitScanForward64(&trailing_zeroes, s);
		s >>= trailing_zeroes;
		for (int_type i = 0; i < 60; i++) {
			int_type val = power_mod(ud(mt), s, p);
			int_type temp = s;
			while (temp != p - 1 and 1 < val and val < p - 1) {
				val = mul_mod(val, val, p);
				temp <<= 1;
			}
			if (val != p - 1 and temp % 2 == 0) return false;
		}
		return true;
	}

	//https://cp-algorithms.com/algebra/factorization.html
	template<typename int_type>
	int_type brent(int_type n, int_type x0 = 2, int_type c = 1) {
		int_type x = x0;
		int_type g = 1;
		int_type q = 1;
		int_type xs, y;

		int_type m = 128;
		int_type l = 1;
		while (g == 1) {
			y = x;
			for (int_type i = 1; i < l; i++)
				x = f(x, c, n);
			int_type k = 0;
			while (k < l && g == 1) {
				xs = x;
				for (int_type i = 0; i < m && i < l - k; i++) {
					x = f(x, c, n);
					q = mul_mod(q, tpp::abs(y,x), n);
				}
				g = std::gcd(q, n);
				k += m;
			}
			l *= 2;
		}
		if (g == n) {
			do {
				xs = f(xs, c, n);
				g = std::gcd(tpp::abs(y, x), n);
			} while (g == 1);
		}
		return g;
	}
		
	template<typename int_type, class mersenne_twister_engine>
	void factor_recursive(int_type n, std::map< int_type, int_type>& previous, mersenne_twister_engine& mt, unsigned long long depth) {
		if (n == 1) return;
		std::uniform_int<int_type> ud(2,std::min(n - 1, int_type(std::numeric_limits<int_fast32_t>::max())));
		auto x0 = ud(mt);
		auto c = ud(mt);
		if (ispmiller(n, mt)) {
			auto it = previous.find(n);
			if (it != previous.end()) //if it's in then either it's fucked or it's ok to do this
				it->second++;
			else previous[n] = 1;
		}
		else {
			int_type result = brent(n, x0, c);
			auto it = previous.find(result);
			if (it != previous.end()) //if it's in then either it's fucked or it's ok to do this
				it->second++;
			else {
				if (result == n || result == 1) {
					if (depth == 0) {// give up
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
	}

	template<typename int_type, class mersenne_twister_engine>
	std::vector<pair<int_type> > factor(int_type n, mersenne_twister_engine& mt) {
		std::vector<pair<int_type> > result;
		std::map<int_type, int_type> factorization_map;
		//factorization_map[1] = 1;
		n = trial_division3(n, factorization_map, int_type(std::numeric_limits<int_fast32_t>::max()));
		if (n > 1) {
			if (ispmiller(n, mt))
				factorization_map[n] = 1;
			else
				factor_recursive(n, factorization_map, mt, 4);
		}
		result.assign(factorization_map.begin(), factorization_map.end());
		return result;
	};

	template<typename int_type>
	std::vector<pair<int_type> > factor(int_type n) {
		std::random_device rd;
		std::mt19937 mt(rd());
		return factor(n, mt);
	};

	template<typename int_type>
	std::tuple<pair<int_type> > euclidian_algorithm(int_type big, int_type small, int_type cutoff) {
		return pair<int_type>(1, 0);
	};

	template<typename int_type>
	std::vector<triple<int_type> > generate_triples(const std::vector<pair<int_type> >& factorization) {
		int_type multiplier = 1;
		std::vector<pair> good_primes;
		for (auto const& it : factorization) {
			if (std::get<0>(it) % 4 == 1) {
				good_primes.push_back(it);
			}
			else {
				multiplier *= ipow(std::get<1>(it), std::get<0>(it));
			}
		}
		std::vector<triple<int_type> > result;
		return result;
	}

	template<typename int_type, class mersenne_twister_engine>
	std::vector<triple<int_type> > generate_triples(int_type n, mersenne_twister_engine& mt) {
		return generate_triples(factor(n, mt));
	}
	template<typename int_type>
	std::vector<triple<int_type> > generate_triples(int_type n) {
		std::random_device rd;
		std::mt19937 mt(rd());
		return generate_triples(factor(n, mt));
	}
}
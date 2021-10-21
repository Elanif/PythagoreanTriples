#pragma once
#include"random.hpp"
#include<map>
#include<numeric>
#include<array>



namespace tpp {
	using Random = effolkronium::random_static;

	typedef std::tuple<int_fast64_t, int_fast64_t> pair;
	typedef std::tuple<int_fast64_t, int_fast64_t, int_fast64_t> triple;

	int64_t ipow(int64_t base, uint8_t exp) { //https://gist.github.com/orlp/3551590
		static const uint8_t highest_bit_set[] = {
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

	long long mult(long long a, long long b, long long mod) {
		long long result = 0;
		while (b) {
			if (b & 1)
				result = (result + a) % mod;
			a = (a + a) % mod;
			b >>= 1;
		}
		return result;
	}

	long long f(long long x, long long c, long long mod) {
		return (mult(x, x, mod) + c) % mod;
	}

	long long brent(long long n, long long x0 = 2, long long c = 1) {
		long long x = x0;
		long long g = 1;
		long long q = 1;
		long long xs, y;

		int m = 64; // 128;
		int l = 1;
		while (g == 1) {
			y = x;
			for (int i = 1; i < l; i++)
				x = f(x, c, n);
			int k = 0;
			while (k < l && g == 1) {
				xs = x;
				for (int i = 0; i < m && i < l - k; i++) {
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

	long long trial_division3(long long n, std::map< int_fast64_t, int_fast64_t>& previous, long long max) {
		auto add = [&](int_fast64_t d) {
			auto it = previous.find(d);
			if (it != previous.end())
				--it->second;
			else
				previous[d] = 1;
		};
		for (int d : {2, 3, 5}) {
			while (n % d == 0) {
				add(d);
				n /= d;
			}
		}
		static std::array<int, 8> increments = { 4, 2, 4, 2, 4, 6, 2, 6 };
		int i = 0;
		for (long long d = 7; d * d <= max; d += increments[i++]) {
			while (n % d == 0) {
				add(d);
				n /= d;
			}
			if (i == 8)
				i = 0;
		}
		return n;
	}

	void factor_recursive(int_fast64_t n, std::map< int_fast64_t, int_fast64_t>& previous, int_fast64_t depth) {
		if (n == 1) return;
		auto x0 = Random::get();
		auto c = Random::get();
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
					factor_recursive(n, previous, depth - 1);
				}
			}
			else {
				factor_recursive(n / result, previous, depth + 1);
				factor_recursive(result, previous, depth + 1);
			}
		}
	}

	std::vector<pair> factor(int_fast64_t n) {
		std::vector<pair> result;
		std::map< int_fast64_t, int_fast64_t> factorization_map;
		//factorization_map[1] = 1;
		n = trial_division3(n, factorization_map, 100000);
		factor_recursive(n, factorization_map, 10);
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

	std::vector<triple> generate_triples(int_fast64_t n) {
		return generate_triples(factor(n));
	}
}
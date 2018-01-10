/*input
10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

/*

************************ Resources ************************

1. Basic: https://www.geeksforgeeks.org/sieve-of-eratosthenes/
2. Segmented Sieve:https://www.geeksforgeeks.org/segmented-sieve/
3. Optimised (Bitwise) Sieve: https://www.quora.com/What-is-bitwise-sieve

*/


std::vector<int> pre() { // calculate all primes upto 10^6
	bool sieve[1000010];
	memset(sieve, true, sizeof(sieve));
	std::vector<int> v;
	sieve[0] = sieve[1] = false;
	for (int i = 2; i <= 1000; i++) {
		if (sieve[i]) {
			for (int j = i+i; j <= 1000000; j += i)
				sieve[j] = false;
		}
	}
	for (int i = 1; i <= 1000000; i++) {
		if (sieve[i])
			v.pb(i);
	}
	return v;
} // O (log(log(n)))

int main() {
	std::vector<int> primes = pre();
	int n; // first n prime numbers (numbers <= 10^6)
	scanf("%d", &n);
	n = min(n, (int) primes.size());
	for (int i = 0; i < n; i++)
		printf("%d ", primes[i]);
	printf("\n");
	return 0;
}
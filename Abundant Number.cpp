/*input
12
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int sum_of_factors(int n) {
	int ans = 0, s = (int) sqrt(n);
	for (int i = 1; i <= s; i++) {
		if (n % i == 0) {
			ans += i;
			if (i != n/i) ans += n/i;
		}
	}
	return ans;
}

bool is_abundant(int n) { // abundant if sum of all proper divisors of number > number
	int s = sum_of_factors(n) - n;
	return (s > n);
}

int main() {
	int n;
	scanf("%d", &n);
	printf("%s\n", is_abundant(n) ? "Yes" : "No");
	return 0;
}
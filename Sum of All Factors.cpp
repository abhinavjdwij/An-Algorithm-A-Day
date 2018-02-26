/*input
120
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

int main() {
	int n;
	scanf("%d", &n);
	printf("Sum of all factors of %d = %d\n", n, sum_of_factors(n));
	return 0;
}
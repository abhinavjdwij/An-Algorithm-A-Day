/*input
2 10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(int n, int k) {
	if (n == 1 || k == 1) return k;
	int ans = INT_MAX;
	for (int i = 1; i <= k; i++) {
		ans = min(ans, max(solve(n-1, i-1), solve(n, k-i)));
	}
	return ans+1;
}

int main() {
	int n, k; // n eggs, k floors (n >= 1 && k >= 1)
	scanf("%d%d", &n, &k);
	printf("%d\n", solve(n, k));
	return 0;
}
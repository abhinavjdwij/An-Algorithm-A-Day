/*input
5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(int n) {
	vector<int> dp (n+1, 0);
	dp[0] = 1;
	for (int i = 1; i < n; i++) {
		for (int j = i; j <= n; j++) {
			dp[j] += dp[j-i];
		}
	}
	return dp[n];
}

int main() {
	int n; // n >= 1 && n <= 10^6
	scanf("%d", &n);
	printf("%d\n", solve(n));
	return 0;
}
/*input
3
4
1 2 3
4 5 1
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

/*

************************ Resources ************************

1. https://www.youtube.com/watch?v=8LusJS5-AGo

*/

int solve(int n, int w, int wt[], int val[]) {
	int dp[110][110];
	memset(dp, 0, sizeof(dp));
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= w; j++) {
			if (j < wt[i]) dp[i][j] = dp[i-1][j];
			else dp[i][j] = max(val[i] + dp[i-1][j-wt[i]], dp[i-1][j]);
		}
	}
	return dp[n][w];
} // O(n*n)

int main() {
	int n, w, wt[110], val[110];
	scanf("%d%d", &n, &w);
	for (int i = 1; i <= n; i++)
		scanf("%d", &val[i]);
	for (int i = 1; i <= n; i++)
		scanf("%d", &wt[i]);
	printf("%d\n", solve(n, w, wt, val));
	return 0;
}
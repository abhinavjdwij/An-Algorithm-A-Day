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

int dp[110][110];

int solve(int w, int p, int n, int wt[], int val[]) {
	if (dp[w][p] != -1) return dp[w][p];
	if (p == n) return 0;
	int ans = (w - wt[p] >= 0) ? val[p] + solve(w - wt[p], p+1, n, wt, val) : -1;
	return dp[w][p] = max(ans, solve(w, p+1, n, wt, val));
} // O(n*n)

int main() {
	int n, w, wt[110], val[110];
	memset(dp, -1, sizeof(dp));
	scanf("%d%d", &n, &w);
	for (int i = 0; i < n; i++)
		scanf("%d", &val[i]);
	for (int i = 0; i < n; i++)
		scanf("%d", &wt[i]);
	printf("%d\n", solve(w, 0, n, wt, val));
	return 0;
}
/*input
4
1 2 3 4 5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;


int n, data[101], dp[101][101]; // ith-matrix has dimension data[i-1] x data[i] for i = 1 to n

int solve(int i, int j) {
	if (dp[i][j] != -1) return dp[i][j];
	if (i == j) return dp[i][j] = 0;
	int k, ans = INT_MAX;
	for (k = i; k < j; k++)
		ans = min(solve(i, k)+solve(k+1, j)+data[i-1]*data[k]*data[j], ans);
	return dp[i][j] = ans;
}

int main() {
	memset(dp, -1, sizeof(dp));
	int n, i;
	scanf("%d", &n);
	for (i = 0; i <= n; i++)
		scanf("%d", &data[i]);
	printf("%d\n", solve(1, n-1));
	return 0;
}

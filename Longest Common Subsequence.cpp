/*input
abbcad
beceafd
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

/*

************************ Resources ************************

1. Video Lecture: https://www.youtube.com/watch?v=NnD96abizww

*/

int dp[110][110];

int lcsq(char a[], char b[], int n, int m) {
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= m; j++) {
			if (a[i-1] == b[j-1])
				dp[i][j] = dp[i-1][j-1] + 1;
			else
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}
	return dp[n][m];
}

int main() {
	int n, m;
	char a[101], b[101];
	scanf("%s%s", a, b);
	n = strlen(a);
	m = strlen(b);
	memset(dp, 0, sizeof(dp));
	printf("%d\n", lcsq(a, b, n, m));
	return 0;
}
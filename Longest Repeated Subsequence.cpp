/*input
aabcbdee
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int lrs(char s[], int n) {
	int dp[110][110];
	memset(dp, 0, sizeof(dp));
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= n; j++) {
			if (s[i-1] == s[j-1] && i != j)
				dp[i][j] = dp[i-1][j-1] + 1;
			else
				dp[i][j] = max(dp[i-1][j], dp[i][j-1]);
		}
	}
	return dp[n][n];
}

int main() {
	int n, m;
	char s[101];
	scanf("%s", s);
	n = strlen(s);
	printf("%d\n", lrs(s, n));
	return 0;
}
/*input
lookingforsomethingoranythinghere
nothinganythingsomething
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int lcs(char a[], int x, char b[], int y) {
	int ans = 0, dp[110][110];
	memset(dp, 0, sizeof(dp));
	for (int i = 1; i <= x; i++) {
		for (int j = 0; j <= y; j++) {
			if (a[i-1] == b[j-1]) {
				dp[i][j] = 1 + dp[i-1][j-1];
				ans = max(ans, dp[i][j]);
			}
			else {
				dp[i][j] = 0;
			}
		}
	}
	return ans;
}

int main() {
	char a[101], b[101];
	scanf("%s\n%s\n", a, b);
	int x = strlen(a), y = strlen(b);
	printf("%d\n", lcs(a, x, b, y));
	return 0;
}
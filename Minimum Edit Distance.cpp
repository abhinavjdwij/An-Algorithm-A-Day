/*input
FOOD
MONEY
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
#define min3(a, b, c) (min(min(a, b), c))
using namespace std;

int solve(string a, string b) {
	int dp[101][101];
	for (int i = 0; i <= a.size(); i++) {
		for (int j = 0; j <= b.size(); j++) {
			if (i == 0) dp[i][j] = j;
			else if (j == 0) dp[i][j] = i;
			else if (a[i-1] != b[j-1]) dp[i][j] = min3(dp[i][j-1], dp[i-1][j], dp[i-1][j-1]) + 1;
			else dp[i][j] = dp[i-1][j-1];
		}
	}
	return dp[a.size()][b.size()];
}

int main() {
	string a, b;
	cin >> a >> b;
	cout << solve(a, b) << endl;
	return 0;
}
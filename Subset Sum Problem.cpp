/*input
5
1 2 3 4 5
8
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

bool subset_sum(int arr[], int n, int s) {
	vector< vector<bool > > dp (n+1, vector<bool> (s+1));
	for (int i = 0; i < n; i++)
		dp[i][0] = true;
	for (int i = 1; i <= s; i++)
		dp[0][i] = false;
	for (int i = 1; i <= n; i++) {
		for (int j = 1; j <= s; j++) {
			if (j < arr[i-1]) dp[i][j] = dp[i-1][j];
			else dp[i][j] = dp[i-1][j] || dp[i-1][j-arr[i-1]];
		}
	}
	return dp[n][s];
}

int main() {
	int n, arr[101], s;
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	scanf("%d", &s);
	cout << (subset_sum(arr, n, s) ? "Yes" : "No");
	return 0;
}
/*input
3 3
0 0 0
0 1 0
0 0 0
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(vector< vector<int> > &arr) {
    if (arr.size() == 0) return 0;
    int m = arr.size(), n = arr[0].size(); 
    vector< vector<int> > dp (m+1, vector<int> (n+1, 0));
    dp[0][1] = 1;
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (arr[i-1][j-1] == 1) dp[i][j] = 0;
            else dp[i][j] = dp[i-1][j] + dp[i][j-1];
        }
    }
    return dp[m][n];
}

int main() {
	int m, n;
	vector< vector<int> > arr; // 1 means obstacle, 0 means empty space
	scanf("%d%d", &m, &n);
	for (int i = 0; i < m; i++) {
		vector<int> data (n);
		for (int i = 0; i < n; i++) {
			scanf("%d", &data[i]);
		}
		arr.pb(data);
	}
	printf("%d\n", solve(arr));
	return 0;
}
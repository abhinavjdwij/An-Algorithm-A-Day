/*input
6 5
0 1 1 0 1
1 1 0 1 0
0 1 1 1 0
1 1 1 1 0
1 1 1 1 1
0 0 0 0 0
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(vector< vector<int> > arr, int n, int m) {
	int ans = 1;
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < m; j++) {
			if (arr[i][j] == 1) {
				arr[i][j] = min(arr[i][j-1], min(arr[i-1][j], arr[i-1][j-1])) + 1;
				ans = max(ans, arr[i][j]);
			}
			else {
				arr[i][j] = 0;
			}
		}
	}
	return ans;
}

int main() {
	int n, m;
	vector< vector<int> > arr;
	scanf("%d%d", &n, &m);
	for (int i = 0; i < n; i++) {
		vector<int> data (m);
		for (int i = 0; i < m; i++) {
			scanf("%d", &data[i]);
		}
		arr.pb(data);
	}
	printf("%d\n", solve(arr, n, m));
	return 0;
}
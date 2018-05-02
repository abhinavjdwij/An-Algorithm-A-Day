/*input
6
10 50 40 30 20 0
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(vector<int> arr, int l, int r) {
	int m = (l+r) / 2, n = arr.size();
	if ((m == 0 || arr[m-1] <= arr[m])
		&& (m == n-1 || arr[m+1] <= arr[m])) {
		return m;
	}
	else if (m > 0 && arr[m-1] > arr[m]) {
		return solve(arr, l, m-1);
	}
	else { return solve(arr, m+1, r); }
} // returns 0-indexed position

int main() {
	int n;
	scanf("%d", &n);
	vector<int> arr (n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &arr[i]);
	}
	printf("%d\n", solve(arr, 0, arr.size()-1));
	return 0;
}
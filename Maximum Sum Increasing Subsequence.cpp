/*input
7
1 101 2 3 100 4 5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(vector<int> arr, int n) {
	vector<int> dp (n);
	for (int i = 0; i < n; i++) {
		dp[i] = arr[i];
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			if (arr[i] > arr[j])
				dp[i] = max(dp[i], arr[i] + dp[j]);
		}
	}
	return *std::max_element(dp.begin(), dp.end());
} // Variation of Longest Increasing Subsequence

int main() {
	int n;
	scanf("%d", &n);
	vector<int> arr (n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	printf("%d\n", solve(arr, n));
	return 0;
}
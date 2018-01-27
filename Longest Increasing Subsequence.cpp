/*input
8
10 22 9 33 21 50 41 60
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int lis(int arr[], int n) {
	int dp[n];
	for (int i = 0; i < n; i++)
		dp[i] = 1;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			if (arr[i] > arr[j])
				dp[i] = max(dp[i], 1 + dp[j]);
		}
	}
	return *std::max_element(dp, dp+n);
} // O(n^2)

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	printf("%d\n", lis(arr, n));
	return 0;
}
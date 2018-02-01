/*input
4
1 2 3 4
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(int arr[], int n) {
	int s = 0, ans = 0;
	for (int i = 0; i < n; i++)
		s += arr[i];
	for (int i = 0; i < n; i++) {
		s -= arr[i];
		ans += s * arr[i];
	}
	return ans;
} // O(n)

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	printf("%d\n", solve(arr, n));
	return 0;
}
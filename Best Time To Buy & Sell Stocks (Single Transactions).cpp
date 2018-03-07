/*input
5
10 2 6 4 8
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(int arr[], int n) {
    if (n < 2) return 0;
    int m = arr[0], ans = arr[1]-arr[0];
    for (int i = 1; i < n; i++) {
        ans = max(ans, arr[i]-m);
        m = min(m, arr[i]);
    }
    return max(ans, 0);
}

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	printf("Max Profit Possible: %d\n", solve(arr, n));
	return 0;
}
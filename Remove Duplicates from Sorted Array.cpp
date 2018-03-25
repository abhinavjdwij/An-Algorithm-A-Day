/*input
10
1 2 2 3 3 3 4 5 5 6
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void solve(int arr[], int &n) {
	int l = 0, r = 0;
	while (r < n) {
		if (arr[r] != arr[l]) {
			l++;
			arr[l] = arr[r];
		}
		r++;
	}
	n = l+1;
}

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	solve(arr, n);
	printf("%d\n", n);
	for (int i = 0; i < n; i++)
		printf("%d ", arr[i]);
	printf("\n");
	return 0;
}
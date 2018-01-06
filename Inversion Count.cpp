/*input
5
2 3 8 6 1
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int merge(int arr[], int l, int m, int r) {
	int i = l, j = m+1, c = 0;
	std::vector<int> temp;
	while (i <= m && j <= r) {
		if (arr[i] < arr[j]) temp.pb(arr[i++]);
		else {
			c += (m-i+1);
			temp.pb(arr[j++]);
		}
	}
	while (i <= m) temp.pb(arr[i++]);
	while (j <= r) temp.pb(arr[j++]);
	for (int i = l, k = 0; i <= r; i++, k++)
		arr[i] = temp[k];
	return c;
}

int invcnt(int arr[], int l, int r) {
	if (l < r) {
		int m = (l+r) / 2;
		return invcnt(arr, l, m) + invcnt(arr, m+1, r) + merge(arr, l, m, r);
	}
	return 0;
}

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	printf("%d\n", invcnt(arr, 0, n-1));
	return 0;
}
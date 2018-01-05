/*input
10
5 2 1 7 3 4 9 8 10 6
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void merge(int arr[], int l, int m, int r) {
	int i = l, j = m+1;
	std::vector<int> temp;
	while (i <= m && j <= r) {
		if (arr[i] < arr[j]) temp.pb(arr[i++]);
		else temp.pb(arr[j++]);
	}
	while (i <= m) temp.pb(arr[i++]);
	while (j <= r) temp.pb(arr[j++]);
	for (int i = l, k = 0; i <= r; i++, k++)
		arr[i] = temp[k];
}

void mergesort(int arr[], int l, int r) {
	if (l < r) {
		int m = (l+r) / 2;
		mergesort(arr, l, m);
		mergesort(arr, m+1, r);
		merge(arr, l, m, r);
	}
}

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	mergesort(arr, 0, n-1);
	for (int i = 0; i < n; i++)
		printf("%d ", arr[i]);
	printf("\n");
	return 0;
}
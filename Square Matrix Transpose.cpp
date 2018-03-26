/*input
4
1 1 1 1
2 2 2 2
3 3 3 3
4 4 4 4
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void transpose(int arr[101][101], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			std::swap(arr[i][j], arr[j][i]);
		}
	}
}

int main() {
	int n, arr[101][101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			scanf("%d", &arr[i][j]);
		}
	}
	transpose(arr, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%d ", arr[i][j]);
		}
		printf("\n");
	}
	return 0;
}
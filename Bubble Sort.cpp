/*input
10
4 6 1 2 7 10 9 5 8 3
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

/*

************************ Resources ************************

1. https://www.youtube.com/watch?v=Jdtq5uKz-w4

*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int bubblesort(int arr[], int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n-i-1; j++) {
			if (arr[j] > arr[j+1])
				std::swap(arr[j], arr[j+1]);
		}
	}
} // O(n^2)

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	bubblesort(arr, n);
	printf("%d\n", n);
	for (int i = 0; i < n; i++)
		printf("%d ", arr[i]);
	printf("\n");
	return 0;
}
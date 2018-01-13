/*input
10
4 6 1 3 7 5 2 10 8 9
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

/*

************************ Resources ************************

1. https://www.youtube.com/watch?v=GUDLRan2DWM

*/


void selectsort(int arr[], int n) {
	for (int i = 0; i < n; i++) {
		int m = i;
		for (int j = i+1; j < n; j++) {
			if (arr[j] < arr[m])
				m = j;
		}
		std::swap(arr[m], arr[i]);
	}
}

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	selectsort(arr, n);
	printf("%d\n", n);
	for (int i = 0; i < n; i++)
		printf("%d ", arr[i]);
	printf("\n");
	return 0;
}
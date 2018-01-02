/*input
10
1 -2 3 4 5 -6 7 -8 9 10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

/*

************************ Resources ************************

1. Basic: https://www.geeksforgeeks.org/largest-sum-contiguous-subarray/
2. Simple Variation: https://www.geeksforgeeks.org/smallest-sum-contiguous-subarray/
3. Advanced 2-D Variation: https://www.geeksforgeeks.org/dynamic-programming-set-27-max-sum-rectangle-in-a-2d-matrix/

*/

int kadane(int arr[], int n) { // Largest Sum Contiguous Subarray
	int m = arr[0], c = arr[0];
	for (int i = 1; i < n; i++) {
		c = max(arr[i], c + arr[i]);
		m = max(c, m);
	}
	return m;
} // O(n)

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	printf("%d\n", kadane(arr, n));
	return 0;
}

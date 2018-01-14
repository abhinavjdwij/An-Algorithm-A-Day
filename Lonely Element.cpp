/*input
11
1 4 2 5 3 4 1 5 2 1 1
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

/*

************************ Resources ************************

1. My Blog: https://medium.com/@dwij28/lonely-element-fd2f813d8614

*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int lonely(int arr[], int n) {
	int ans = 0;
	for (int i = 0; i < n; i++)
		ans ^= arr[i];
	return ans;
} // O(n)

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	printf("%d\n", lonely(arr, n));
	return 0;
}
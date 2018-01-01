/*input
10
1 2 3 4 4 4 5 6 6 7
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

/*

************************ Resources ************************

1. Basic: https://www.geeksforgeeks.org/binary-search/
2. Variations: https://www.geeksforgeeks.org/the-ubiquitous-binary-search-set-1/
3. Interesting Alternate Approach: http://codeforces.com/blog/entry/9901
4. Advanced Tutorial: https://www.topcoder.com/community/data-science/data-science-tutorials/binary-search/

*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int search(int arr[], int l, int r, int x) { // 0-indexed position of x in [l, r] or -1 if it does not occur in [l, r]
	int m;
	while (l <= r) {
		m = (l+r) / 2;
		if (arr[m] == x) return m;
		else if (arr[m] < x) l = m+1;
		else r = m-1;
	}
	return -1;
} // O(log(n))

int search_recursive(int arr[], int l, int r, int x) { // 0-indexed position of x in [l, r] or -1 if it does not occur in [l, r]
	if (l > r) return -1;
	int m = (l+r) / 2;
	if (arr[m] == x) return m;
	else if (arr[m] < x) return search_recursive(arr, m+1, r, x);
	else return search_recursive(arr, l, m-1, x);
} // O(log(n))

int lowerbound(int arr[], int l, int r, int x) { // 0-indexed first occurrence position of x in [l, r] or (r+1) if it does not occur
	int m, p = r+1;
	while (l <= r) {
		m = (l+r) / 2;
		if (arr[m] == x) {
			p = m;
			r = m-1;
		}
		else if (arr[m] < x) l = m+1;
		else r = m-1;
	}
	return p;
} // O(log(n))

int upperbound(int arr[], int l, int r, int x) { // 0-indexed last occurrence position of x in [l, r] or (r+1) if it does not occur
	int m, p = r+1;
	while (l <= r) {
		m = (l+r) / 2;
		if (arr[m] == x) {
			p = m;
			l = m+1;
		}
		else if (arr[m] < x) l = m+1;
		else r = m-1;
	}
	return p;
} // O(log(n))

int count_occurrences(int arr[], int n, int x) { // Number of occurrences of x in arr
	int lb = lowerbound(arr, 0, n-1, x);
	int ub = upperbound(arr, 0, n-1, x);
	if (lb == n) return 0;
	else return ub - lb + 1;
} // O(log(n))

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	printf("%d\n", search(arr, 0, n-1, 5));
	printf("%d\n", search_recursive(arr, 0, n-1, 5));
	printf("%d\n", search(arr, 0, n-1, 10));
	printf("%d\n", search_recursive(arr, 0, n-1, 10));
	printf("%d\n", lowerbound(arr, 0, n-1, 6));
	printf("%d\n", upperbound(arr, 0, n-1, 6));
	printf("%d\n", count_occurrences(arr, n, 6));
	printf("%d\n", lowerbound(arr, 0, n-1, 4));
	printf("%d\n", upperbound(arr, 0, n-1, 4));
	printf("%d\n", count_occurrences(arr, n, 4));
	printf("%d\n", lowerbound(arr, 0, n-1, 5));
	printf("%d\n", upperbound(arr, 0, n-1, 5));
	printf("%d\n", count_occurrences(arr, n, 5));
	printf("%d\n", lowerbound(arr, 0, n-1, 9));
	printf("%d\n", upperbound(arr, 0, n-1, 9));
	printf("%d\n", count_occurrences(arr, n, 9));
	return 0;
}
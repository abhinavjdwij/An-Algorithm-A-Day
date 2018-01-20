/*input
9
1 2 3 4 5 6 7 8 9
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int main() {
	int n, arr[110], l, r;
	scanf("%d", &n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &arr[i]);
	}
	l = 0;
	r = n-1;
	while (l < r) {
		while (arr[l] % 2 == 0 && l < r) // increment l if even numbers found
			l++;
		while (arr[r] % 2 && l < r) // decrement r if odd numbers found
			r--;
		if (l < r) {
			std::swap(arr[l], arr[r]);
			l++;
			r--;
		}
	}
	for (int i = 0; i < n; i++)
		printf("%d ", arr[i]);
	printf("\n");
	return 0;
}
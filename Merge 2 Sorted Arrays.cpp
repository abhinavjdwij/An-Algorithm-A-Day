/*input
5
1 3 4 6 10
5
2 5 7 8 9
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

/*

************************ Resources ************************

1. Used in Merge Sort: https://www.youtube.com/watch?v=KF2j-9iSf4Q

*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int merge(int a, int x[], int b, int y[], int ans[]) {
	int i = 0, j = 0, k = 0;
	while (i < a && j < b) {
		if (x[i] <= y[j]) ans[k++] = x[i++];
		else ans[k++] = y[j++];
	}
	while (i < a) ans[k++] = x[i++];
	while (j < b) ans[k++] = y[j++];
	return k;
} // O(a+b)

int main() {
	int a, b, x[101], y[101], n, ans[202];
	scanf("%d", &a);
	for (int i = 0; i < a; i++)
		scanf("%d", &x[i]);
	scanf("%d", &b);
	for (int i = 0; i < b; i++)
		scanf("%d", &y[i]);
	n = merge(a, x, b, y, ans);
	printf("%d\n", n);
	for (int i = 0; i < n; i++)
		printf("%d ", ans[i]);
	return 0;
}
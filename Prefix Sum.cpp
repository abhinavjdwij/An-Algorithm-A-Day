/*input
10
1 2 3 4 5 6 7 8 9 10
5
1 5
2 5
3 10
1 10
4 7
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int main() {
	int n, q, l, r, arr[11], data[11+1];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	data[0] = 0;
	for (int i = 1; i <= n; i++)
		data[i] = data[i-1] + arr[i-1];
	scanf("%d", &q);
	while (q--) {
		scanf("%d%d", &l, &r);
		printf("%d\n", data[r] - data[l-1]);
	}
	return 0;
}
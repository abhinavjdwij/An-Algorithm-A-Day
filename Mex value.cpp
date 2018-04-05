/*input
5
1 0 3 5 4
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int mex(int data[], int n) {
	int ans = 0, p = 0;
	sort(data, data+n);
	while (p < n) {
		if (data[p] != ans) return ans;
		ans++;
		p++;
	}
	return ans;
} // minimum excluded positive integer

int main() {
	int n, data[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &data[i]);
	}
	printf("%d\n", mex(data, n));
	return 0;
}
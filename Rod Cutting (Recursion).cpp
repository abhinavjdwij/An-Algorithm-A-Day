/*input
8
1 5 8 9 10 17 17 20
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(int n, vector<int> price) {
	if (n <= 0) return 0;
	int ans = INT_MIN;
	for (int i = 0; i < n; i++) {
		ans = max(ans, price[i] + solve(n-i-1, price));
	}
	return ans;
}

int main() {
	int n;
	scanf("%d", &n);
	vector<int> price (n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &price[i]);
	}
	printf("%d\n", solve(n, price));
	return 0;
}
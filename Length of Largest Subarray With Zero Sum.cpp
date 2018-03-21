/*input
8
15 -2 2 -8 1 7 10 23
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(vector<int> arr) {
	map<int, int> data;
	data[0] = 0;
	int s = 0, ans = 0;
	for (int i = 0; i < arr.size(); i++) {
		s += arr[i];
		if (data.count(s)) ans = max(ans, (i+1) - data[s]);
		else data[s] = (i+1);
	}
	return ans;
}

int main() {
	int n;
	vector<int> arr;
	scanf("%d", &n);
	for (int i = 0; i < n; i++) {
		int x;
		scanf("%d", &x);
		arr.pb(x);
	}
	printf("%d\n", solve(arr));
	return 0;
}
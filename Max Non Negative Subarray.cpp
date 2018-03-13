/*input
6
1 2 5 -7 2 3
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

vector<int> solve(vector<int> &arr) {
    long long cs = 0, s = 0;
    int cl = 0, cr = -1, l = 0, r = -1;
    for (int i = 0; i < arr.size(); i++) {
        if (arr[i] < 0) { cl = i+1; cs = 0; }
        else { cs += arr[i]; cr = i; }
        if (cs > s || (cs == s && cr-cl > r-l)) { l = cl; r = cr; s = cs; }
    }
    vector<int> ans;
    if (r >= 0 && r < arr.size()) {
        for (int i = l; i <= r; i++)
            ans.push_back(arr[i]);
    }
    return ans;
}

int main() {
	int n;
	scanf("%d", &n);
	vector<int> arr;
	for (int i = 0; i < n; i++) {
		int x;
		scanf("%d", &x);
		arr.pb(x);
	}
	vector<int> ans = solve(arr);
	printf("%ld\n", ans.size());
	for (int i : ans) printf("%d ", i);
	printf("\n");
	return 0;
}
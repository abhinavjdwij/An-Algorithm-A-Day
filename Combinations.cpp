/*input
4 2
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void solve(vector< vector<int> > &ans, vector<int> &arr, int n, int k) {
    if (arr.size() > k) return;
    if (arr.size() == k) ans.push_back(arr);
    int p = 1;
    if (arr.size() > 0) p = arr[arr.size()-1] + 1;
    for (int i = p; i <= n; i++) {
        arr.push_back(i);
        solve(ans, arr, n, k);
        arr.pop_back();
    }
}

vector< vector<int> > combinations(int n, int k) {
    vector< vector<int> > ans;
    if (n < k) return ans;
    vector<int> arr;
    solve(ans, arr, n, k);
    return ans;
} // return all possible combinations of k numbers out of 1 2 3 ... n

int main() {
	int n, k;
	scanf("%d%d", &n, &k);
	vector< vector<int> > ans = combinations(n, k);
	for (auto v : ans) {
		for (int i : v) {
			printf("%d ", i);
		}
		printf("\n");
	}
	return 0;
}
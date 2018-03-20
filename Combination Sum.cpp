/*input
4
2 3 6 7
7
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void solve(set <vector<int> > &ans, vector<int> arr, vector<int> &data, int p, int k) {
    if (k == 0) ans.insert(data);
    else {
        for (int i = p; i < arr.size() && k >= arr[i]; i++) {
            data.pb(arr[i]);
            solve(ans, arr, data, i, k-arr[i]);
            data.pop_back();
        }
    }
}

set< vector<int> > combinationSum(vector<int> &arr, int k) {
    std::sort(arr.begin(), arr.end());
    set< vector<int> > ans;
    vector<int> data;
    solve(ans, arr, data, 0, k);
    return ans;
}

int main() {
    int n, k;
    vector<int> arr;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) {
        int x;
        scanf("%d", &x);
        arr.pb(x);
    }
    scanf("%d", &k);
    set< vector<int> > ans = combinationSum(arr, k);
    for (auto s : ans) {
        for (int i : s) {
            printf("%d ", i);
        }
        printf("\n");
    }
    return 0;
}
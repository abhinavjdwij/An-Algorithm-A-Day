/*input
4
2 7 11 15
9
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

pair<int, int> solve(const vector<int> &arr, int s) { // returns 0-based index of solution
    pair<int, int> ans;
    if (arr.size() < 2) return ans;
    map<int, int> data;
    for (int i = 0; i < arr.size(); i++) {
        if (data.find(s - arr[i]) != data.end()) {
            ans = mp(data[s - arr[i]], i);
            return ans;
        }
        else if (data.find(arr[i]) == data.end()) {
            data[arr[i]] = i;
        }
    }
    return ans;
}

int main() {
	int n, s;
	vector<int> arr;
	pair<int, int> ans;
	scanf("%d", &n);
	for (int i = 0; i < n; i++) {
		int x;
		scanf("%d", &x);
		arr.pb(x);
	}
	scanf("%d", &s);
	ans = solve(arr, s);
	printf("%d + %d = %d\n", arr[ans.first], arr[ans.second], s);
	return 0;
}
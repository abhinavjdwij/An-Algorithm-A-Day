/*input
6
2 1 5 6 2 3
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(vector<int> &arr) {
    stack<int> s;
    int ans = 0, p = 0;
    while (p < arr.size()) {
        if (s.empty() || arr[s.top()] <= arr[p]) {
            s.push(p++);
        }
        else {
            int t = s.top();
            s.pop();
            ans = max(ans, arr[t]* (s.empty() ? p : p - s.top()-1));
        }
    }
    while (!s.empty()) {
        int t = s.top();
        s.pop();
        ans = max(ans, arr[t]* (s.empty() ? p : p - s.top()-1));
    }
    return ans;
}

int main() {
	int n;
	scanf("%d", &n);
	vector<int> arr (n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &arr[i]);
	}
	printf("%d\n", solve(arr));
	return 0;
}
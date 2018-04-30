/*input
5
1 2 + 5 *
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(int l, int r, string s) {
    if (s == "+") return (r+l);
    else if (s == "-") return (r-l);
    else if (s == "*") return (r*l);
    else return (r/l);
}

int evaluate(vector<string> &arr) {
    stack<int> data;
    for (string s : arr) {
        if (s == "+" || s == "-" || s == "*" || s == "/") {
            int l = data.top();
            data.pop();
            int r = data.top();
            data.pop();
            data.push(solve(l, r, s));
        }
        else {
            data.push(atoi(s.c_str()));
        }
    }
    return data.top();
}

int main() {
	int n;
	cin >> n;
	vector<string> arr(n);
	for (int i = 0; i < n; i++) {
		cin >> arr[i];
	}
	cout << evaluate(arr) << endl;
	return 0;
}
/*input
1
51
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

string solve(int n) {
	string ans = "";
	while (n) {
		int r = n % 26;
		if (r == 0) {
			ans += 'Z';
			n = n/26 - 1;
		}
		else {
			ans += 'A' + (r-1);
			n /= 26;
		}
	}
	std::reverse(ans.begin(), ans.end());
	return ans;
}

int main() {
	int t, n;
	cin >> t;
	while (t--) {
		cin >> n;
		cout << solve(n) << endl;
	}
	return 0;
}
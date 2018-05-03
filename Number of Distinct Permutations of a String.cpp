/*input
abbacbdb
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

vector<int> fac;

void pre() {
	fac.pb(1);
	for (int i = 1; i <= 10; i++) {
		fac.pb(i*fac[i-1]);
	}
}

int solve(string s) {
	map<char, int> data;
	for (char c : s) {
		data[c]++;
	}
	int k = 1;
	for (auto i : data) {
		k *= fac[i.second];
	}
	int ans = fac[s.size()] / k;
	return ans;
}


int main() {
	pre();
	string s;
	cin >> s;
	cout << solve(s) << endl;
	return 0;
}
/*input
AbhiNAVJha
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

string lowercase(string s) {
	string ans = "";
	for (char i : s) {
		if (i >= 'A' && i <= 'Z') {
			i += 32;
		}
		ans += i;
	}
	return ans;
}

int main() {
	string s;
	cin >> s;
	cout << lowercase(s) << endl;
	return 0;
}
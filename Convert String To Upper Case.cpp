/*input
AbhiNAVJha
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

string uppercase(string s) {
	string ans = "";
	for (char i : s) {
		if (i >= 'a' && i <= 'z') {
			i -= 32;
		}
		ans += i;
	}
	return ans;
}

int main() {
	string s;
	cin >> s;
	cout << uppercase(s) << endl;
	return 0;
}
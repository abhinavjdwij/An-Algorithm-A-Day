/*input
ThIs iS A SENtenCe.
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

string sentencecase(string s) {
	string ans = "";
	for (char i : s) {
		if (i >= 'A' && i <= 'Z') {
			i += 32;
		}
		ans += i;
	}
	if (ans.size()) {
		if (ans[0] >= 'a' && ans[0] <= 'z') {
			ans[0] -= 32;
		}
	} 
	return ans;
}

int main() {
	string s;
	getline(cin, s);
	cout << sentencecase(s) << endl;
	return 0;
}
/*input
AbhinavJha
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

string toggleCase(string s) {
	for (int i = 0; i < s.size(); i++) {
		s[i] ^= 32;
	}
	return s;
}

int main() {
	string s;
	cin >> s;
	cout << toggleCase(s) << endl;
	return 0;
}
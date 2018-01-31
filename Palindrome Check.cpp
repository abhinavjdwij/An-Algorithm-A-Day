/*input
abccba
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

bool pal(string s) {
	int n = s.size();
	for (int i = 0; i <= n/2; i++) {
		if (s[i] != s[n-i-1])
			return false;
	}
	return true;
}

int main() {
	string s;
	cin >> s;
	cout << (pal(s) ? "Palindrome" : "Non-Palindrome") << endl;
	return 0;
}
/*input
AD
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(string s) {
    int ans = 0;
    for (std::string::iterator i = s.begin(); i != s.end(); i++)
        ans = ans*26 + (*i - 64);
    return ans;
}

int main() {
	string s;
	cin >> s;
	cout << solve(s);
	return 0;
}
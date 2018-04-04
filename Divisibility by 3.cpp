/*input
123456789
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

bool divby3(string s) {
	int k = 0;
	for (char i : s) {
		k += (int) (i - '0');
	}
	return (k % 3 == 0);
}

int main() {
	string s; // large number as string
	cin >> s;
	cout << ((divby3(s)) ? "Yes" : "No") << endl;
	return 0;
}
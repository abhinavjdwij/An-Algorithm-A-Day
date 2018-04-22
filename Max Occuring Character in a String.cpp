/*input
samplestring
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

char max_occur(string s) {
	char ans = (char) 0;
	int m = 0;
	map<char, int> data;
	for (char i : s) {
		data[i]++;
		if (data[i] > m) {
			m = data[i];
			ans = i;
		}
	}
	return ans;
}

int main() {
	string s;
	cin >> s;
	cout << max_occur(s) << endl;
	return 0;
}
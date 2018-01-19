/*input
abcbadecdab
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void counter(string s) {
	std::map<char, int> data;
	for (int i = 0; i < s.size(); i++)
		data[s[i]]++;
	for (auto i : data)
		cout << i.first << " occurs " << i.second << " times." << endl;
}

int main() {
	string s;
	cin >> s;
	counter(s);
	return 0;
}
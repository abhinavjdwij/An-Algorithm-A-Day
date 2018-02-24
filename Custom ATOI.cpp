/*input
-123
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int my_atoi(string str) {
    int x = 0, flag = 1;
    for (int i = 0; i < str.size(); i++) {
        if (str[i] == '-') flag = -1;
        else if (str[i] >= '0' && str[i] <= '9') x = x*10 + (str[i]-'0');
        else return -1;
    }
    return x*flag;
}

int main() {
	string s;
	int n;
	cin >> s;
	n = my_atoi(s);
	cout << typeid(s).name() << ": " << s << endl;
	cout << typeid(n).name() << ": " << n << endl;
	return 0;
}
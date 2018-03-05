/*input
10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

string bin(int n) {
    if (n == 0) return "0";
    string ans = "";
    while (n) {
        ans += (n % 2) + '0';
        n /= 2;
    }
    std::reverse(ans.begin(), ans.end());
    return ans;
}

int main() {
	int n;
	cin >> n;
	cout << bin(n);
	return 0;
}
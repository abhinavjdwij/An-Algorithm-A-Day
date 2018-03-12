/*input
5
3 30 34 5 9
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int compare(string a, string b) { return a+b > b+a; }

string largestNumber(vector<int> &arr) {
    std::vector<string> data;
    for (int i = 0; i < arr.size(); i++)
        data.push_back(to_string(arr[i]));
    std::sort(data.begin(), data.end(), compare);
    string ans = "";
    for (int i = 0; i < data.size(); i++)
        ans += data[i];
    if (ans.size() && ans[0] == '0') return "0"; 
    return ans;
}


int main() {
	int n;
	std::vector<int> data;
	cin >> n;
	for (int i = 0; i < n; i++) {
		int x;
		cin >> x;
		data.pb(x);
	}
	cout << largestNumber(data);
	return 0;
}
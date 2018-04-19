/*input
CXIV
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int val(char c) {
    switch(c) {
        case 'I': return 1;
        case 'V': return 5;
        case 'X': return 10;
        case 'L': return 50;
        case 'C': return 100;
        case 'D': return 500;
        case 'M': return 1000;
        default: return 0;
    }
    return -1;
}

int convert(string s) {
    vector<int> arr;
    for (char i : s) {
        arr.push_back(val(i));
    }
    int ans = 0;
    for (int i = 0; i < arr.size(); i++) {
        if (i+1 < arr.size())
            if (arr[i] >= arr[i+1]) {
                ans += arr[i];
            }
            else {
                ans += arr[i+1] - arr[i];
                i++;
            }
        else {
            ans += arr[i];
        }
    }
    return ans;
}

int main() {
	string s;
	cin >> s;
	cout << convert(s) << endl;
	return 0;
}
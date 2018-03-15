/*input
81
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int sqroot(int n) {
    int l = 0, r = 46340, m, ans;
    while (l <= r) {
        m = (l+r) / 2;
        if (m*m == n) {
            ans = m;
            break;
        }
        else if (m*m < n) {
            l = m+1;
            ans = m;
        }
        else r = m-1;
    }
    return ans;
}

int main() {
	int n;
	scanf("%d", &n);
	printf("%d\n", sqroot(n));
	return 0;
}
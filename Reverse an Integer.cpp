/*input
245
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int reverse(int n) {
    bool flag = (n < 0);
    n = abs(n);
    int ans = 0;
    while (n) {
        ans = ans*10 + (n % 10);
        n /= 10;
    }
    if (flag) ans *= -1;
    return ans;
}

int main() {
	int n;
	scanf("%d", &n);
	printf("%d\n", reverse(n));
	return 0;
}
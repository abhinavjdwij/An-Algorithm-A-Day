/*input
10 3
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(int a, int b) {
    bool flag = ((a >= 0 && b >= 0) || (a < 0 && b < 0));
    long long x = abs((long long) a);
    long long y = abs((long long) b);
    long long z = 1;
    while (y < x) {
        y = y << 1;
        z = z << 1;
    }
    long long ans = 0;
    while (x >= abs((long long) b)) {
        while (y <= x) {
            x -= y;
            ans += z;
        }
        y = y >> 1;
        z = z >> 1;
    }
    if (!flag) ans *= -1;
    return ((ans > INT_MAX || ans < INT_MIN) ? INT_MAX : ans);
} // returns floor of result and checks overflows

int main() {
	int a, b;
	scanf("%d%d", &a, &b);
	printf("%d\n", solve(a, b));
	return 0;
}
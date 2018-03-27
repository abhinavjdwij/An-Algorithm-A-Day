/*input
3 5
1 10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void solve(int a, int b, int c, int d, int &x, int &y) {
	y = (b*d) / __gcd(b, d);
	x = a*(y/b) + c*(y/d);
} // sum (a/b) + (c/d);

int main() {
	int a, b, c, d, x, y;
	scanf("%d%d", &a, &b);
	scanf("%d%d", &c, &d);
	solve(a, b, c, d, x, y);
	printf("%d/%d\n", x, y);
	return 0;
}
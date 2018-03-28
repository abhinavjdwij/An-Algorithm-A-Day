/*input
6 10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void solve(int &a, int &b) {
	int k = __gcd(a, b);
	a /= k;
	b /= k;
}

int main() {
	int a, b; // a/b
	scanf("%d%d", &a, &b);
	solve(a, b);
	printf("%d/%d\n", a, b);
	return 0;
}
/*input
18
24
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int gcd(int a, int b) {
	while (a) {
		int t = b % a;
		b = a;
		a = t;
	}
	return b;
}

int main() {
	int a, b;
	scanf("%d%d", &a, &b);
	printf("%d\n", gcd(a, b));
	return 0;
}
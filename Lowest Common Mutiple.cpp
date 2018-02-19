/*input
12
18
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

int lcm(int a, int b) {
	return a*b/gcd(a, b);
}

int main() {
	int a, b;
	scanf("%d%d", &a, &b);
	printf("%d\n", lcm(a, b));
	return 0;
}
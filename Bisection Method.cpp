/*input
1 100 130
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

float f (float x) { return x*x - 81; }

int main() {
	int n;
	float a, b, c, root;
	scanf("%f%f%d", &a, &b, &n); // lower limit, upper limit and number of iterations 
	if (! f(a)) printf("Root = %f\n", a);
	else if (! f(b)) printf("Root = %f\n", b);
	else {
		for (int i = 1; i <= n; i++) {
			c = (a+b)/2;
			if (! f(c)) {
				printf("Root = %f\n", c);
				exit(0);
			}
			else if (f(a)*f(c) < 0) { b = c; }
			else if (f(b)*f(c) < 0) { a = c; }
		}
		printf("Approximate Root = %f\n", c);
	}
	return 0;
}
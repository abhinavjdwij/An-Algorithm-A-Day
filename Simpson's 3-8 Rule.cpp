/*input
4 5.2 6
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

float f (float x) { return log(x) ; }

int main() {
	int i, n;
	float a, b, h, s, multiple3 = 0, others = 0;
	scanf("%f%f%d", &a, &b, &n); // lower limit, upper limit and no. of intervals
	h = (b-a)/n;
	if (n % 3 != 0) printf("Simpson's 3/8 Rule is not applicable.\n");
	else {
		for (i = 1; i <= n-1; i++) {
			if (i % 3 == 0) multiple3 += f(a + i*h);
			else others += f(a + i*h);
		}
		s = (3 * h / 8) * (f(a) + f(b) + 3*others + 2*multiple3);
		printf("Integral = %f\n", s);
	}
	return 0;
}
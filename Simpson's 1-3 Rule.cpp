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
	float a, b, h, s, odd = 0.0, even = 0.0;
	scanf("%f%f%d", &a, &b, &n); // lower limit, upper limit and no. of intervals
	h = (b-a)/n;
	if (n % 2 != 0) printf("Simpson's 1/3 Rule is not applicable.\n");
	else {
		for (i = 1; i <= n-1; i++) {
			if (i % 2 == 0) even += f(a + i*h);
			else odd += f(a + i*h);
		}
		s = (h / 3) * (f(a) + f(b) + 2*even + 4*odd);
		printf("Integral = %f\n", s);
	}
	return 0;
}
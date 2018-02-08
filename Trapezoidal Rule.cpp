/*input
0 1
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

float f (float x) { return 4*x - 3*x*x; } // function

int main() {
	int n;
	float a, b, h, s;
	scanf("%f%f%d", &a, &b, &n); // lower limit, upper limit and no. of intervals
	h = (b-a)/n;
	s = f(a) + f(b);
	for (int i = 1; i < n; i++)
		s += 2*f(a + i*h);
	printf("Integral = %f\n", h*s/2);
	return 0;
}
/*input
0 0 0.05 0.15
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

float f(float x, float y) { return x*x + y*y; }

int main() {
	float x, y, h, eval;
	scanf("%f%f%f%f", &x, &y, &h, &eval); // x0, y0, step size & point of evaluation
	while (x < eval) {
		y = y + h*f(x, y);
		x += h;
	}
	printf("Result = %.4f\n", y);
	return 0;
}
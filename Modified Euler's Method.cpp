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
	float x[100], y[100], h, eval, diff, temp;
	scanf("%f%f%f%f", &x[0], &y[0], &h, &eval); // x0, y0, step size & point of evaluation
	int i = 0;
	for (i = 1; x[i-1] < eval; i++) {
		x[i] = x[i-1] + h;
		y[i] = y[i-1] + h*f(x[i-1], y[i-1]);
		diff = 100.0; // huge value;
		temp = y[i];
		while (diff > 0.0001) {
			temp = y[i-1] + h * (f(x[i-1], y[i-1]) + f(x[i], y[i])) / 2;
			diff = temp - y[i];
			y[i] = temp;
		}
	}
	printf("Result = %f\n", y[i-1]);
	return 0;
}
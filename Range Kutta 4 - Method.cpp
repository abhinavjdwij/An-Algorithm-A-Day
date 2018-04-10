/*input
0 1 0.1 0.4
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

float f (float x, float y) { return x-y; }

int main() {
	float x, y, h, eval, k1, k2, k3, k4;
	scanf("%f%f%f%f", &x, &y, &h, &eval); // x0, y0, h and point of evaluation
	while (x < eval) {
		k1 = h*f(x, y);
		k2 = h*f(x + (h/2.0), y + (k1/2.0));
		k3 = h*f(x + (h/2.0), y + (k2/2.0));
		k4 = h*f(x + h, y + k3);
		y += (k1 + 2*k2 + 2*k3 + k4) / 6;
		x += h;
	}
	printf("Answer = %f\n", y);
	return 0;
}
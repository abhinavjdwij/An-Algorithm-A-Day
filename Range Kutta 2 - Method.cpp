/*input
0 1 0.1 0.4
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

float f (float x, float y) {  return x-y; }

int main() {
	float x, y, h, eval, k1, k2;
	scanf("%f%f%f%f", &x, &y, &h, &eval); // x0, y0, h and point of evaluation
	while (x < eval) {
		k1 = h*f(x, y);
		k2 = h*f(x+h, y+k1);
		y += 0.5 * (k1 + k2);
		x += h;
	}
	printf("Answer = %f\n", y);
	return 0;
}
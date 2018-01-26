/*input
1.6 10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

float f (float x) { return 3*x - cos(x) - 1; }

float df(float x) { return 3 + sin(x); } 

int main() {
	int i, n;
	float x;
	scanf("%d%f", &n, &x); // number of iterations and initial value of root
	for (i = 1; i <= n; i++) {
		x = x - (f(x)/df(x));
	}
	printf("Value = %f\n", x);
	return 0;
}
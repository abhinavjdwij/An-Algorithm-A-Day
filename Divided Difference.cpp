/*input
4
0 2
1 3
2 12
5 147
4
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int main() {
	float D[10], x[10], x1, res, v = 1.0;
	int i, j, n;
	scanf("%d", &n); // number of terms
	for (i = 0; i < n; i++)
		scanf("%f%f", &x[i], &D[i]); // x & y values
	scanf("%f", &x1); // point of evaluation
	res = D[0];
	for (i = 1; i < n; i++) {
		for (j = 0; j < n-i; j++) {
			D[i] = (D[j] - D[j+1]) / (x[i+j] - x[j+1]);
		}
		v *= (x1 - x[i-1]);
		res += v*D[0];
	}
	printf("Result = %f\n", res);
	return 0;
}
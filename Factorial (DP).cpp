/*input
5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int fac[11];

int pre() {
	fac[0] = fac[1] = 1;
	for (int i = 2; i <= 10; i++)
		fac[i] = fac[i-1]*i;
}

int main() {
	pre();
	int n;
	scanf("%d", &n);
	printf("fac(%d) = %d\n", n, fac[n]);
	return 0;
}
/*input
10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int fib[51];

int pre() {
	fib[1] = 0;
	fib[2] = 1;
	for (int i = 3; i <= 40; i++)
		fib[i] = fib[i-1] + fib[i-2];
}

int main() {
	pre();
	int n;
	scanf("%d", &n);
	printf("fib(%d) = %d\n", n, fib[n]);
	return 0;
}
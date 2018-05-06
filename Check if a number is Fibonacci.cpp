/*input
21
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int isPerfectSq(int n) {
	int s = (int) sqrt(n);
	return (n == s*s);
}

int isFibo(int n) {
	int x = 5*n*n + 4, y = 5*n*n - 4;
	return (isPerfectSq(x) || isPerfectSq(y));
}

int main() {
	int n;
	scanf("%d", &n);
	printf("%s\n", isFibo(n) ? "Yes" : "No");
	return 0;
}
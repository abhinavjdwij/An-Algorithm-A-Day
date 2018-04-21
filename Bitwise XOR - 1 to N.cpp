/*input
5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int xor1toN(int n) {
	int k = n % 4;
	if (k == 0) return n;
	else if (k == 1) return 1;
	else if (k == 2) return n+1;
	else return 0;
}

int main() {
	int n;
	scanf("%d", &n);
	printf("%d\n", xor1toN(n));
	return 0;
}
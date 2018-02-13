/*input
9
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int phi(int n) {
	int ans = n;
	for (int i = 2; i*i <= n; i++) {
		if (n % i == 0) {
			while (n % i == 0) {
				n /= i;
			}
			ans -= ans/i;
		}
	}
	if (n > 1) ans -= ans/n;
	return ans;
}

int main() {
	int n;
	scanf("%d", &n);
	printf("%d\n", phi(n));
	return 0;
}
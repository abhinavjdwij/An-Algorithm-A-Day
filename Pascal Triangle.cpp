/*input
5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void pascal(int n) {
	for (int i = 1; i <= n; i++) {
		int x = 1;
		for (int j = 1; j <= i; j++) {
			printf("%d ", x);
			x = x*(i-j)/j;
		}
		printf("\n");
	}
}

int main() {
	int n;
	scanf("%d", &n);
	pascal(n);
	return 0;
}
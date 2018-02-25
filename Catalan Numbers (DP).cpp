/*input
10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int catalan[25];

void pre() {
	memset(catalan, 0, sizeof(catalan));
	catalan[1] = 1; catalan[2] = 1;
	for (int i = 3; i <= 20; i++) {
		for (int j = 1; j <= i; j++)
			catalan[i] += catalan[j]*catalan[i-j];
	}
}

int main() {
	pre();
	int n;
	scanf("%d", &n);
	printf("%d\n", catalan[n]);
	return 0;
}
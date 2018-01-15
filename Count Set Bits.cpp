/*input
23
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int csb(int n) { // count number of '1's in binary representation of n
	int c = 0;
	while (n) {
		n &= (n-1);
		c++;
	}
	return c;
} // O(log(n))

int main() {
	int n;
	scanf("%d", &n);
	printf("%d\n", csb(n));
	return 0;
}
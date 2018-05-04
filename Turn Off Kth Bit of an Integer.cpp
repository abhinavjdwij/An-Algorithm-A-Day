/*input
15 3
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int turnOffKthBit(int n, int k) {
	return (n & ~(1 << (k - 1)));
} // k >= 1

int main() {
	int n, k;
	scanf("%d%d", &n, &k);
	printf("%d\n", turnOffKthBit(n, k));
	return 0;
}
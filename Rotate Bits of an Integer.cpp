/*input
16 2
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int left_rotation(int n, int k) {
	return ((n << k) | (n >> (32-k)));
}

int right_rotation(int n, int k) {
	return ((n >> k) | (n << (32-k)));
}

int main() {
	int n, k; // rotate n by k bits
	scanf("%d%d", &n, &k);
	printf("Left Rotation of %d by %d bits = %d\n", n, k, left_rotation(n, k));
	printf("Right Rotation of %d by %d bits = %d\n", n, k, right_rotation(n, k));
	return 0;
}
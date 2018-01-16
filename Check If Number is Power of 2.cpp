/*input
16
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

bool check(int n) { return (n && ! (n & (n-1))); } // O(1)

int main() {
	int n;
	scanf("%d", &n);
	printf("%s\n", check(n) ? "Yes" : "No");
	return 0;
}
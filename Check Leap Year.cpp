/*input
2018
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

bool check(int n) { return (n % 4 == 0 && n % 100 != 0 || n % 400 == 0); }

int main() {
	int n;
	scanf("%d", &n);
	printf("%s\n", check(n) ? "Leap" : "Non-Leap Year");
	return 0;
}
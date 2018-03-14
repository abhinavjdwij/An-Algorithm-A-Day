/*input
abcde
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int main() {
	int n;
	char s[101];
	scanf("%s", s);
	n = strlen(s);
	for (int i = 0; i < n/2; i++)
		swap(s[i], s[n-i-1]);
	printf("%s\n", s);
	return 0;
}
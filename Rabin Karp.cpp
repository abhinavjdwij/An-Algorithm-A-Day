/*input
LookingForSomethingHere
Something
LookingForSomethingHere
Nothing
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

/*

************************ Resources ************************

1. https://www.youtube.com/watch?v=H4VrKHVG5qI

*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

bool check(char *a, char *b, int p) {
	for (int i = 0; b[i] != '\0'; i++)
		if (a[p+i] != b[i])
			return false;
	return true;
}

bool rabin_karp(char *text, char *s) {
	int hash = 0, c = 0, x = strlen(s), y = strlen(text);
	if (x > y) return false;
	for (int i = 0; i < x; i++) {
		hash += s[i];
		c += text[i];
	}
	if (hash == c && check(text, s, 0)) return true;
	for (int i = x; i < y; i++) {
		c += text[i] - text[i-x];
		if (hash == c && check(text, s, i-x+1)) return true;
	}
	return false;
}

int main() {
	char text[110], s[110];
	scanf("%s%s", text, s);
	printf("%s\n", rabin_karp(text, s) ? "Yes" : "No");
	scanf("%s%s", text, s);
	printf("%s\n", rabin_karp(text, s) ? "Yes" : "No");
	return 0;
}
/*input
aabbcccdefff
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

char* encode(char *src) {
    int p = 0, c = 1;
    char *s;
    s = (char *) malloc(sizeof(char)*110);
    for (int i = 1; src[i] != '\0'; i++) {
        if (src[i] == src[i-1]) {
            c += 1;
        }
        else {
            s[p++] = src[i-1];
            s[p++] = '0' + c;
            c = 1;
        }
    }
    s[p++] = src[strlen(src)-1];
    s[p++] = '0' + c;
    s[p] = '\0';
    return s;
}

int main() {
    char s[101];
    scanf("%s", s);
    printf("%s\n", encode(s));
    return 0;
}
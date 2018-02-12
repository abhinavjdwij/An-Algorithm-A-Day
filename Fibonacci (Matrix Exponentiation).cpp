/*input
10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void multiply(long long int f[2][2], long long int m[2][2]) {
    long long int a = (f[0][0] * m[0][0]) + (f[0][1] * m[1][0]);
    long long int b = (f[0][0] * m[0][1]) + (f[0][1] * m[1][1]);
    long long int c = (f[1][0] * m[0][0]) + (f[1][1] * m[1][0]);
    long long int d = (f[1][0] * m[0][1]) + (f[1][1] * m[1][1]);
    f[0][0] = a; f[0][1] = b;
    f[1][0] = c; f[1][1] = d;
}

void modExpo(long long int f[2][2], long long int n) {
    if (n == 0 || n == 1) return;
    long long int m[2][2] = {{1, 1}, {1, 0}};
    modExpo(f, n/2);
    multiply(f, f);
    if (n % 2 != 0) multiply(f, m);
}

long long int fib(long long int n) {
    long long int f[2][2] = {{1, 1}, {1, 0}};
    if (n == 0) return 0;
    modExpo(f, n-1);
    return f[0][0];
}

int main() {
    long long int n;
    scanf("%lld", &n);
    printf("%lld", fib(n));
    return 0;
}
/*input
3
1 4 5
2
2 3
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

double solve(const vector<int> &a, int la, int ra, const vector<int> &b, int lb, int rb, int k) {
    if (la > ra) return b[k+lb];
    if (lb > rb) return a[k+la];
    if (k < 1) return min(a[k+la], b[k+lb]);
    int i = (la+ra) / 2, j = (lb+rb) / 2;
    int x = a[i], y = b[j];
    if ((i-la)+(j-lb) < k) {
        if (x > y)
            return solve(a, la, ra, b, j+1, rb, k-(j-lb)-1);
        else
            return solve(a, i+1, ra, b, lb, rb, k-(i-la)-1);
    }
    else {
        if (x > y)
            return solve(a, la, i-1, b, lb, rb, k);
        else
            return solve(a, la, ra, b, lb, j-1, k);
    }
}

double median(vector<int> &a, vector<int> &b) {
    int x = a.size(), y = b.size();
    int k = x+y;
    if (k % 2 == 0)
        return (solve(a, 0, x-1, b, 0, y-1, k/2) + solve(a, 0, x-1, b, 0, y-1, k/2 - 1)) / 2.0;
    else
        return solve(a, 0, x-1, b, 0, y-1, k/2);
} // O (log(m+n))

int main() {
	int m, n, x;
	vector<int> a;
	vector<int> b;
	scanf("%d", &m);
	for (int i = 0; i < m; i++) {
		scanf("%d", &x);
		a.pb(x);
	}
	scanf("%d", &n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &x);
		b.pb(x);
	}
	printf("%lf\n", median(a, b));
	return 0;
}
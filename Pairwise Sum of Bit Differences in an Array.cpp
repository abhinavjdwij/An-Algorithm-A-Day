/*input
3
1 3 5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(vector<int> &arr) {
    long long ans = 0, mod = 1000000007, n = arr.size();
    for (int i = 0; i < 32; i++) {
        long long c = 0;
        for (int p : arr) {
            if (p & (1 << i)) {
                c++;
            }
        }
        long long k = (c*(n-c)*2);
        ans = (ans + k) % mod;
    }
    return ((int) ans);
}

int main() {
	int n;
	scanf("%d", &n);
	vector<int> arr (n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &arr[i]);
	}
	printf("%d\n", solve(arr));
	return 0;
}
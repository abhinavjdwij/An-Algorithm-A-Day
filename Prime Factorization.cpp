/*input
90
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

std::vector<int> prfac(int n) {
	std::vector<int> v;
	if (n < 2) return v;
	int s = (int) sqrt(n);
	if (n % 2 == 0) {
		v.pb(2);
		while (n % 2 == 0)
			n /= 2;
	}
	for (int i = 3; i <= s; i += 2) {
		if (n % i == 0) {
			v.pb(i);
			while (n % i == 0)
				n /= i;
		}
	}
	if (n > 1) v.pb(n);
	return v;
} // O(sqrt(n))

int main() {
	int n;
	scanf("%d", &n);
	std::vector<int> ans = prfac(n);
	printf("Prime factors of %d: ", n);
	for (int i = 0; i < ans.size(); i++)
		printf("%d ", ans[i]);
	printf("\n");
	return 0;
}
/*input
120
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

std::vector<int> factors(int n) {
	std::vector<int> data;
	int s = (int) sqrt(n);
	for (int i = 1; i <= s; i++) {
		if (n % i == 0) {
			data.pb(i);
			if (i != n/i) data.pb(n/i);
		}
	}
	std::sort(data.begin(), data.end());
	return data;
}

int main() {
	int n;
	scanf("%d", &n);
	std::vector<int> data = factors(n);
	printf("Factors of %d:\n", n);
	for (auto i : data)
		printf("%d ", i);
	printf("\n");
	return 0;
}
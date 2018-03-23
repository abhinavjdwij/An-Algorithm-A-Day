/*input
10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

vector<int> solve(int n) {
	vector<int> data(n, 0);
	for (int i = 1; i <= n; i++) {
		data[i] = (3*i*i - i) / 2;
	}
	return data;
}

int main() {
	int n;
	scanf("%d", &n);
	vector<int> data = solve(n);
	printf("The first %d Pentagonal Numbers are: ", n);
	for (int i = 1; i <= n; i++)
		printf("%d ", data[i]);
	printf("\n");
	return 0;
}

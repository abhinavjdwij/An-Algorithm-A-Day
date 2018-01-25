/*input
5
1 2 3 4 5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

std::vector< std::vector<int> > gensub(std::vector<int> arr, int n) {
	std::vector< std::vector<int> > ans;
	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			std::vector<int> temp;
			for (int k = i; k <= j; k++) {
				temp.pb(arr[k]);
			}
			ans.pb(temp);
		}
	}
	return ans;
}

int main() {
	int n;
	std::vector<int> arr;
	std::vector< std::vector<int> > ans;
	scanf("%d", &n);
	for (int i = 0, x = 0; scanf("%d", &x) && i < n; i++)
		arr.pb(x);
	ans = gensub(arr, n);
	for (int i = 0; i < ans.size(); i++) {
		for (int j = 0; j < ans[i].size(); j++)
			printf("%d ", ans[i][j]);
		printf("\n");
	}
	return 0;
}
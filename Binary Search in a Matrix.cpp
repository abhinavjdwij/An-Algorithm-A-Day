/*input
4 3
1 3 5 7
10 11 16 20
23 30 34 50
34
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

bool bs(std::vector<int> data, int x) {
    int l = 0, r = data.size()-1, m;
    while (l <= r) {
        m = (l+r) / 2;
        if (data[m] < x) l = m+1;
        else if (data[m] > x) r = m-1;
        else return true;
    }
    return false;
}

bool bs_matrix(std::vector<vector<int> > &arr, int x) {
    int i = 0, j = 0, n = arr.size(), m = arr[0].size();
    bool ans = false;
    while (i < n) {
        if (arr[i][m-1] == x) return true;
        else if (arr[i][m-1] < x) i++;
        else return bs(arr[i], x);
    }
    return ans;
}

int main() {
	int n, m, x;
	std::vector<vector<int> > arr;
	scanf("%d%d", &n, &m);
	for (int i = 0; i < n; i++) {
		std::vector<int> temp;
		for (int j = 0; j < m; j++) {
			scanf("%d", &x);
			temp.pb(x);
		}
		arr.pb(temp);
	}
	scanf("%d", &x);
	printf("%s\n", bs_matrix(arr, x) ? "Yes" : "No");
	return 0;
}
/*input
4
2
3 4
6 5 7
4 1 8 3
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int solve(vector<vector<int> > &arr) {
    if (arr.size() == 0) return 0;
    for (int i = arr.size()-2; i >= 0; i--) {
        for (int j = 0; j < arr[i].size(); j++) {
            arr[i][j] += min(arr[i+1][j], arr[i+1][j+1]);
        }
    }
    return arr[0][0];
} // movement allowed only to adjacent numbers on the row below

int main() {
	int n;
	vector< vector<int> > arr;
	scanf("%d", &n);
	for (int i = 0; i < n; i++) {
		vector<int> data (i+1);
		for (int j = 0; j <= i; j++) {
			scanf("%d", &data[j]);
		}
		arr.pb(data);
	}
	printf("%d\n", solve(arr));
	return 0;
}
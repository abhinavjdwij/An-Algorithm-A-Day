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

void solve(int arr[], int n) {
	map< int, pair<int, int> > data;
	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			int s = arr[i] + arr[j];
			if (data.find(s) != data.end()) {
				cout << arr[data[s].first] << " + " << arr[data[s].second] << " = " << arr[i] << " + " << arr[j];
				return;
			}
			else { data[s] = mp(i, j); }
		}
	}
	cout << "No Equal Pairs Exist" << endl;
}

int main() {
	int n, arr[101];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	solve(arr, n);
	return 0;
}
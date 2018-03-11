/*input
3
1 2 3
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void solve(int p, vector<int> data, int arr[], int n) {
    if (p == n) {
    	for (int i : data)
    		printf("%d ", i);
    	printf("\n");
    }
    else {
        solve(p+1, data, arr, n);
        data.push_back(arr[p]);
        solve(p+1, data, arr, n);
        data.pop_back();
    }
}

void allsubsets(int arr[], int n) {
    vector<int> data;
    solve(0, data, arr, n);
}

int main() {
	int n, arr[11];
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &arr[i]);
	allsubsets(arr, n);
	return 0;
}
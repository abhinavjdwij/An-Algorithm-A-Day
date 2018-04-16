/*input
6
4 1 1 2 1 3
1
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void remove_all(vector<int> &arr, int x) {
    int l = 0, r = 0;
    while (r < arr.size()) {
        if (arr[r] != x) arr[l++] = arr[r];
        r++;
    }
    arr.resize(l);
}

int main() {
	int n, x;
	scanf("%d", &n);
	vector<int> arr (n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &arr[i]);
	}
	scanf("%d", &x);
	remove_all(arr, x);
	n = arr.size();
	printf("%d\n", n);
	for (int i : arr) {
		printf("%d ", i);
	}
	printf("\n");
	return 0;
}
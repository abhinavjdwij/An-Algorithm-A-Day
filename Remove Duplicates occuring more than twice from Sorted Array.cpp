/*input
11
1 1 1 1 2 2 3 3 3 4 5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void solve(vector<int> &arr) {
    if (arr.size() <= 2) return;
    int i = 0, j = 1;
    bool flag = true;
    while (j < arr.size()) {
        if (arr[i] != arr[j]) {
            i++;
            arr[i] = arr[j];
            flag = true;
        }
        else if (flag) {
            i++;
            arr[i] = arr[j];
            flag = false;
        }
        j++;
    }
    arr.resize(i+1);
}

int main() {
	int n;
	scanf("%d", &n);
	vector<int> arr (n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &arr[i]);
	}
	solve(arr);
	n = arr.size();
	printf("%d\n", n);
	for (int i : arr) {
		printf("%d ", i);
	}
	printf("\n");
	return 0;
}
/*input
10 3
10 10 9 8 7 6 5 4 3 2
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

vector<int> sliding_window_max(vector<int> &arr, int k) {
    deque<int> data;
    int p = 0;
    vector<int> ans;
    k = min(k, (int) arr.size());
    while (p < k) {
        while ((!data.empty()) && (arr[p] >= arr[data.back()]))
            data.pop_back();
        data.push_back(p);
        p++;
    }
    while (p < arr.size()) {
        ans.push_back(arr[data.front()]);
        while ((!data.empty()) && (data.front() <= p-k))
            data.pop_front();
        while ((!data.empty()) && (arr[p] >= arr[data.back()]))
            data.pop_back();
        data.push_back(p);
        p++;
    }
    ans.push_back(arr[data.front()]);
    return ans;
}

int main() {
	int n, k;
	scanf("%d%d", &n, &k);
	vector<int> arr (n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &arr[i]);
	}
	vector<int> ans = sliding_window_max(arr, k);
	for (int i : ans) {
		printf("%d ", i);
	}
	printf("\n");
	return 0;
}
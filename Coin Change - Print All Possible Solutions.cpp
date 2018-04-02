/*input
3
1 3 5
10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void solve(int w, int arr[], int n, int p, vector<int> data) {
    if (p >= n) return;
    if (w < 0) return;
    if (w == 0) {
        for (int i : data)
            printf("%d ", i);
        printf("\n");
        return;
    }
    if (w >= arr[p]) {
        data.push_back(arr[p]);
        solve(w - arr[p], arr, n, p, data);
        data.pop_back();
    }
    solve(w, arr, n, p+1, data);
}

int main() {
    int n, w, arr[101];
    scanf("%d", &n);
    for (int i = 0; i < n; i++)
        scanf("%d", &arr[i]);
    scanf("%d", &w);
    vector<int> data;
    solve(w, arr, n, 0, data);
    return 0;
}
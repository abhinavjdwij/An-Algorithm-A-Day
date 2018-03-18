/*input
4
4 2 1 3
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

vector<int> wave(vector<int> arr) {
    std::sort(arr.begin(), arr.end());
    for (int i = 0; i < arr.size()-1; i += 2) {
        std::swap(arr[i], arr[i+1]);
    }
    return arr;
} // sort such that a1 >= a2 <= a3 >= a4 <= a5 .... and so on

int main() {
	int n;
	vector<int> arr;
	scanf("%d", &n);
	for (int i = 0; i < n; i++) {
		int x;
		scanf("%d", &x);
		arr.pb(x);
	}
	arr = wave(arr);
	printf("%d\n", n);
	for (int i : arr)
		printf("%d ", i);
	printf("\n");
	return 0;
}
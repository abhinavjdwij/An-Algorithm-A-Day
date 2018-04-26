/*input
8
2 -1 4 3 5 -1 3 2
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

vector<int> longest_increasing_subseq(vector<int> arr) {
	int n = arr.size();
	vector<int> lis (n, 1);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			if (arr[i] > arr[j]) {
				lis[i] = max(lis[i], 1 + lis[j]);
			}
		}
	}
	return lis;
}

vector<int> longest_decreasing_subseq(vector<int> arr) {
	std::reverse(arr.begin(), arr.end());
	vector<int> lds = longest_increasing_subseq(arr);
	std::reverse(arr.begin(), arr.end());
	std::reverse(lds.begin(), lds.end());
	return lds;
}

int longest_bitonic_subseq(vector<int> arr) {
	int n = arr.size();
	vector<int> lis = longest_increasing_subseq(arr);
	vector<int> lds = longest_decreasing_subseq(arr);
	vector<int> lbs (n);
	for (int i = 0; i < n; i++) {
		lbs[i] = lis[i] + lds[i] - 1;
	}
	return *std::max_element(lbs.begin(), lbs.end());
}

int main() {
	int n;
	scanf("%d", &n);
	vector<int> arr (n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &arr[i]);
	}
	printf("%d\n", longest_bitonic_subseq(arr));
	return 0;
}
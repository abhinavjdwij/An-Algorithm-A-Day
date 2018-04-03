/*input
3
3 4 5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

bool solve(int arr[], int n) {
	int k = arr[0];
	for (int i = 1; i < n; i++)
		k ^= arr[i];
	return (k != 0);
}

void predict(int data[], int n) {
	bool ans = solve(data, n);
	printf("%s will win.\n", (ans) ? "Player - 1" : "Player - 2");
} // 2 player nim game predictor where Player - 1 starts first

int main() {
	int n, data[101]; // number of piles and coins in each pile
	scanf("%d", &n);
	for (int i = 0; i < n; i++)
		scanf("%d", &data[i]);
	predict(data, n);
	return 0;
}
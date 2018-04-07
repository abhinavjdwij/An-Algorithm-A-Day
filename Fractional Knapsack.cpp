/*input
3 50
60 10
100 20
120 30
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

bool compare(pair<int, int> a, pair<int, int> b) {
	float x = (float) (a.first / a.second);
	float y = (float) (b.first / b.second);
	return (x > y);
}

float solve(vector< pair<int, int> > arr, int w) {
	std::sort(arr.begin(), arr.end(), compare);
	int c = 0;
	float ans = 0.0;
	for (int i = 0; i < arr.size(); i++) {
		if (c + arr[i].second <= w) {
			c += arr[i].second;
			ans += arr[i].first;
		}
		else {
			float r = (float) (w - c);
			ans += arr[i].first * (r / arr[i].second);
			break;
		}
	}
	return ans;
}

int main() {
	int n, w; // max allowed weight
	vector< pair<int, int> > arr; // first element in pair is value and second is weight
	scanf("%d%d", &n, &w);
	for (int i = 0; i < n; i++) {
		int val, wt;
		scanf("%d%d", &val, &wt);
		arr.push_back(mp(val, wt));
	}
	printf("%f\n", solve(arr, w));
	return 0;
}
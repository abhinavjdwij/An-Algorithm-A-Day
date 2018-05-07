/*input
3
1 2 4
2
3 5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void mergeInPlace(vector<int> &a, vector<int> &b) {
    int i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] > b[j]) {
            a.insert(a.begin()+i, b[j]);
            j++;
        }
        i++;
    }
    while (j != b.size()) {
        a.push_back(b[j]);
        j++;
    }
}

int main() {
	int n, m;
	scanf("%d", &n);
	vector<int> a (n);
	for (int i = 0; i < n; i++) {
		scanf("%d", &a[i]);
	}
	scanf("%d", &m);
	vector<int> b (m);
	for (int i = 0; i < m; i++) {
		scanf("%d", &b[i]);
	}
	mergeInPlace(a, b);
	printf("Inplace merged vector:\n");
	for (int i : a) {
		printf("%d ", i);
	}
	printf("\n");
	return 0;
}
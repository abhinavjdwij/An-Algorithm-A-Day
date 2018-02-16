/*input
5
1 2 100
2 1 19
3 2 27
4 1 25
5 3 15
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

typedef struct job_ {
	int id, deadline, profit;
}job;

int compare(const void *a, const void *b) { // descending
	return (((job *) a)->profit < ((job *) b)->profit);
}

void schedule(job data[], int n) {
	int i, j, check[101], ans = 0;
	memset(check, 0, sizeof(check));
	for (i = 0; i < n; i++) {
		for (j = min(data[i].deadline, n)-1; j >= 0; j--) {
			if (!check[j]) {
				check[j] = data[i].id;
				ans += data[i].profit;
				break;
			}
		}
	}
	printf("Sequence of Jobs: ");
	for (i = 0; i < n; i++)
		if (check[i])
			printf("%d ", check[i]);
	printf("\n");
	printf("Max profit: %d\n", ans);
}

int main() {
	job data[101];
	int n;
	scanf("%d", &n); // number of jobs
	for (int i = 0; i < n; i++)
		scanf("%d%d%d", &data[i].id, &data[i].deadline, &data[i].profit); // jobs in the order (id deadline profit)
	qsort(data, n, sizeof(job), compare);
	schedule(data, n);
	return 0;
}
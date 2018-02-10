/*input
0 1 0
2 0 3
4 0 10
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int main() {
	int matrix[3][3], i, j, x = 0, y = 0;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			scanf("%d", &matrix[i][j]);
	for (i = 0; i < 3; i++) {
			x += matrix[0][i % 3] * matrix[1][(i+1) % 3] * matrix[2][(i+2) % 3];
			y += matrix[2][i % 3] * matrix[1][(i+1) % 3] * matrix[0][(i+2) % 3];
	}
	printf("Value = %d\n", x - y);
	return 0;
}
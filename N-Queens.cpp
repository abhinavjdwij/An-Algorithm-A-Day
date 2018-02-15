/*input
4
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int n;
int board[101][101]; // 1 means Queen is placed, else not

int check(int r, int c) {
	int i, j;
	for (i = 0; i < c; i++) // same row
		if (board[r][i])
			return 0;
	for (i = r, j = c; i >= 0 && j >= 0; i--, j--) // upper left diagonal
		if (board[i][j])
			return 0;
	for (i = r, j = c; i < n && j >= 0; i++, j--) // lower left diagonal
		if (board[i][j])
			return 0;
	return 1;
}

int solve(int col) { // places 1 queen in each column first
	int i, j;
	if (col >= n) return 1;
	for (i = 0; i < n; i++) {
		if (check(i, col)) {
			board[i][col] = 1;
			if (solve(col+1))
				return 1;
			board[i][col] = 0;
		}
	}
	return 0;
}

int main() {
	scanf("%d", &n);
	solve(0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			printf("%d ", board[i][j]);
		}
		printf("\n");
	}
	return 0;
}
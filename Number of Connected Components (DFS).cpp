/*input
10 6
1 2
3 1
4 5
6 7
8 5
9 6
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int n, m; // number of vertices, number of edges
std::vector<int> adj[101]; // adjacency list (1-indexed)
bool vis[101]; // visited check (1-indexed)

void dfs(int x) {
	vis[x] = true;
	for (int i = 0; i < adj[x].size(); i++) {
		if (! vis[adj[x][i]]) {
			dfs(adj[x][i]);
		}
	}
}

int solve() {
	int c = 0;
	memset(vis, false, sizeof(vis));
	for (int i = 1; i <= n; i++) {
		if (! vis[i]) {
			dfs(i);
			c++;
		}
	}
	return c;
}

int main() {
	scanf("%d%d", &n, &m);
	for (int i = 1; i <= m; i++) {
		int u, v;
		scanf("%d%d", &u, &v);
		adj[u].pb(v);
		adj[v].pb(u);
	}
	printf("Number of connected components = %d\n", solve());
	return 0;
}
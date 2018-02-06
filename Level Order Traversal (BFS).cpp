/*input
10 9
1 2
1 3
2 4
2 5
2 6
3 7
7 8
8 9
8 10
1
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int n, m; // number of nodes, number of edges
std::vector<int> adj[101]; // adjacency list
int level[101]; // levels w.r.t root (also used as check for visited nodes)

void bfs(int root) {
	memset(level, -1, sizeof(level));
	std::queue<int> q;
	q.push(root);
	level[root] = 0;
	while (! q.empty()) {
		int x = q.front();
		q.pop();
		for (int i = 0; i < adj[x].size(); i++) {
			if (level[adj[x][i]] == -1) {
				level[adj[x][i]] = 1 + level[x];
				q.push(adj[x][i]);
			}
		}
	}
}

int main() {
	scanf("%d%d", &n, &m);
	for (int i = 1; i <= m; i++) {
		int u, v;
		scanf("%d%d", &u, &v);
		adj[u].pb(v);
		adj[v].pb(u);
	}
	int root; // root of the graph;
	scanf("%d", &root);
	bfs(root);
	for (int i = 1; i <= n; i++) {
		printf("Node %d is at level: %d\n", i, level[i]);
	}
	return 0;
}
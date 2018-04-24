/*input
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

class Graph {
	int n;
	std::vector< pair<int, int> > *adj;
public:
	Graph(int n);
	void add(int u, int v, int w);
	void dijkstra(int src);
};

Graph::Graph(int n) {
	this->n = n;
	adj = new std::vector< pair<int, int> > [n];
}

void Graph::add(int u, int v, int w) {
	adj[v].pb(mp(u, w));
	adj[u].pb(mp(v, w));
}

void Graph::dijkstra(int src) {
	priority_queue< pair<int, int>, std::vector< pair <int, int> >, greater< pair<int, int> > > pq;
	std::vector<int> d(n, INT_MAX);
	pq.push(mp(0, src));
	d[src] = 0;
	while (! pq.empty()) {
		int u = pq.top().second;
		pq.pop();
		for (int i = 0; i < adj[u].size(); i++) {
			if (d[adj[u][i].first] > d[u] + adj[u][i].second) {
				d[adj[u][i].first] = d[u] + adj[u][i].second;
				pq.push(mp(d[adj[u][i].first], adj[u][i].first));
			}
		}
	}
	printf("Vertex  |  Distance From Source\n_______________________________\n");
	for (int i = 0; i < n; i++)
		printf("  %d\t    |    \t%d\n", i, d[i]);
}

int main() {
	int n = 9;
    Graph g(n);
    g.add(0, 1, 4);
    g.add(0, 7, 8);
    g.add(1, 2, 8);
    g.add(1, 7, 11);
    g.add(2, 3, 7);
    g.add(2, 8, 2);
    g.add(2, 5, 4);
    g.add(3, 4, 9);
    g.add(3, 5, 14);
    g.add(4, 5, 10);
    g.add(5, 6, 2);
    g.add(6, 7, 1);
    g.add(6, 8, 6);
    g.add(7, 8, 7);
    g.dijkstra(0);
	return 0;
}
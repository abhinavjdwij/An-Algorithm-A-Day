/*input
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

typedef struct TreeNode {
	int val;
	TreeNode *left;
	TreeNode *right;
	TreeNode(int x) : val(x), left(NULL), right(NULL) {}
}TreeNode;

TreeNode* solve(const vector<int> &arr, int l, int r) {
    if (l > r) return NULL;
    int m = (l+r) / 2;
    TreeNode* p = new TreeNode(arr[m]);
    p->left = solve(arr, l, m-1);
    p->right = solve(arr, m+1, r);
    return p;
}

TreeNode* sortedArrayToBST(const vector<int> &arr) {
    TreeNode* ans = solve(arr, 0, arr.size()-1);
    return ans;
}

int main() {
	return 0;
}
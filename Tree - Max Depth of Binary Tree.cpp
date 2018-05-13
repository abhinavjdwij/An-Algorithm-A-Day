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
 
int Solution::maxDepth(TreeNode* root) {
    if (root == NULL) return 0;
    else if (root->left == NULL && root->right == NULL) return 1;
    else if (root->left == NULL) return 1 + Solution::maxDepth(root->right);
    else if (root->right == NULL) return 1 + Solution::maxDepth(root->left);
    else return 1 + max(Solution::maxDepth(root->left), Solution::maxDepth(root->right));
}

int main() {
	return 0;
}
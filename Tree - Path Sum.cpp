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

int Solution::hasPathSum(TreeNode* root, int k) {
    if (root == NULL) return 0;
    else if (root->val == k && root->left == NULL && root->right == NULL) {
        return 1;
    }
    return max(Solution::hasPathSum(root->left, k - root->val), Solution::hasPathSum(root->right, k - root->val));
}

int main() {
	return 0;
}
/*input
{([])}
(()}]
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

bool isPair(char a, char b) {
	if ((a == '{' && b == '}') || (a == '}' && b == '{')) return true;
	if ((a == '(' && b == ')') || (a == ')' && b == '(')) return true;
	if ((a == '[' && b == ']') || (a == ']' && b == '[')) return true;
	return false;
}

bool balanced(char s[]) {
	char a, b;
	std::stack<char> data;
	for (int i = 0; s[i] != '\0'; i++) {
		if (!data.empty() && isPair(s[i], data.top())) data.pop();
		else data.push(s[i]);
	}
	while (data.size() >= 2) {
		a = data.top(); data.pop();
		b = data.top(); data.pop();
		if (!isPair(a, b)) {
			data.push(a);
			data.push(b);
			break;
		}
	}
	return data.empty();
}

int main() {
	char x[110], y[110];
	scanf("%s", x);
	scanf("%s", y);
	printf("%s : %s\n", x, balanced(x) ? "balanced" : "not balanced");
	printf("%s : %s\n", y, balanced(y) ? "balanced" : "not balanced");
	return 0;
}
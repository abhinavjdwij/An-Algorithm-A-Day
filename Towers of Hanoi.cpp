/*input
4
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

void towers(int n, char from, char to, char aux) { // takes in arguments: no. of disks, initial peg, destination peg(final peg), and auxiliary peg(helper peg)
    if (n == 1) {
        printf("Move disk 1 from %c to %c.\n", from, to);
        return;
    }
    towers(n-1, from, aux, to);
    printf("Move disk %d from %c to %c.\n", n, from, to);
    towers(n-1, aux, to, from);
}

int main() {
    int n;
    long long int moves;
    printf("Enter number of disks: ");
    scanf("%d", &n);
    moves = (long long int) pow(2, n) - 1;
    printf("\nThe number of moves is %lld.\n", moves);
    printf("\nThe sequence of moves is:\n\n");
    towers(n, 'A', 'C', 'B');
    return 0;
}
/*input
10
4 7 1 3 2 6 9 10 8 5
*/

/*~ @author = dwij28 (Abhinav Jha) ~*/

#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
using namespace std;

int partition(int *data, int begin, int end) {
    int i, pivot = data[end],  index = begin; // index will store the final position pivot element goes to, i.e the point of partition
    for (i = begin; i <= end-1; i++) {
        if (data[i] <= pivot) {
            std::swap(data[i], data[index]);
            index++;
        }
    }
    std::swap(data[index], data[end]);
    return index;
}

void quicksort(int *data, int begin, int end) {
    if (begin < end) { // begin == end is the case of single element and begin > end is an invalid case that occurs when one side of the partition is empty
        int index = partition(data, begin, end);
        quicksort(data, begin, index-1);
        quicksort(data, index+1, end);
    }
}

int main() {
    int n, data[101], i;
    scanf("%d", &n);
    for (i = 0; i < n; i++)
        scanf("%d", &data[i]);
    quicksort(data, 0, n-1);
    printf("%d\n", n);
    for (i = 0; i < n; i++)
        printf("%d ", data[i]);
    printf("\n");
    return 0;
}
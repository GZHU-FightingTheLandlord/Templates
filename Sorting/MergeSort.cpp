namespace MergeSort{
#define MAXN 100005
    void msort(int *arr, int l, int r) {
        if (l >= r) return;
        int mid = (l + r) >> 1;
        msort(arr, l, mid);
        msort(arr, mid + 1, r);
        if (arr[mid] <= arr[mid + 1]) return;

        int i = l, j = mid + 1, k = l;
        static int tmp[MAXN];
        while (i <= mid && j <= r) {
            if (arr[i] <= arr[j]) tmp[k++] = arr[i++];
            else tmp[k++] = arr[j++];
        }
        while (i <= mid) tmp[k++] = arr[i++];
        while (j <= r) tmp[k++] = arr[j++];
        while (l <= r) arr[l] = tmp[l], l++;
    }
#undef MAXN
}
using MergeSort::msort;

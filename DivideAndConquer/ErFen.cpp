int search(int arr[], int l, int r, int tar)
{
    int mid;
    while (l < r) {
        mid = (l + r) >> 1;
        if (arr[mid] >= tar) r = mid;
        else l = mid + 1;
    }
    return l;
}

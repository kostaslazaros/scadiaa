import hashlib


def list2md5(lst):
    sorted_lst = [i for i in lst]
    sorted_lst.sort()
    lst2str = "".join(sorted_lst)
    print(lst2str)
    hex =  hashlib.md5(lst2str.encode()).hexdigest()[:15]
    return f'A{hex}'


if __name__ == "__main__":
    lst = ['AD00103', 'AD00102', 'AD00101']
    print(list2md5(lst))
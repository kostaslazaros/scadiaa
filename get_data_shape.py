import gzip


def get_data_shape(file_path):
    with gzip.open(file_path, 'rb') as f:
        first_line = f.readline()
        cells = first_line.count(b'\t')
        genes = sum(1 for _ in f)
    # print("Number of genes:", num_rows)
    # print("Number of cells:", num_tabs)
    return cells, genes

def get_shapes(file_list):
    shape_list = ['id,cells,genes']
    for file_name in file_list:
        cl,gn = get_data_shape(f"{file_name}_expr.txt.gz")
        # shape_list.append(get_data_shape(file_name))
        shape_list.append(f"{file_name},{cl},{gn}")
    with open("shape_list.csv", "w") as f:
        f.write("\n".join(shape_list))


def get_file_list():
    fname_list = []
    with open('scREAD_description.csv', 'r') as f:
        
        for i, line in enumerate(f.readlines()):
            if i == 0:
                continue
            fname, *_ = line.split(',')
            fname_list.append(fname)
    return fname_list


if __name__ == "__main__":
    lst = get_file_list()
    # get_data_shape('AD00110_expr.txt.gz')
    get_shapes(lst)
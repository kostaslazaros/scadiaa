import scanpy as sc
import hdf5plugin
import argparse


def check_h5ad(filename):
    anndata = sc.read_h5ad(filename)
    # print(anndata)
    print(f'matrix shape (cells, genes) = {anndata.X.shape}')
    # print(anndata.var)
    print(f'{anndata.obs.head(1)}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Checking h5ad file')
    parser.add_argument('--file', help='h5ad file name')
    args = parser.parse_args()
    print(args)
    # check_h5ad('./data/h5ad/AD00501.h5ad')
    check_h5ad(args.file)
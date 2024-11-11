import os
import sys
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import scipy.io as sio
import scanpy.external as sce
import matplotlib.pyplot as plt

sc.settings.verbosity = 1

# Adjust Scanpy figure defaults
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,
    facecolor = 'white', figsize=(6,6), format='png')


def sparse_std(X,blocksize=100):
    """
    Return the std deviation of the columns of a sparse matrix.
    
    inputs
        X: sparse matrix

        blocksize: number of columns to calculate in parallel. Larger
        blocksize will increase memory usage as this function converts 
        each block to a dense matrix and uses the numpy .std() function
    
    outputs
        X_std: numpy vector with std deviation of each column of X
    
    """
    n,m = X.shape
    blocksize = 100
    k = int(m/blocksize)
    X_std = []
    for i in range(k):
        X_std += list(X[:,blocksize*i:blocksize*(i+1)].todense().std(0).A1)
    X_std += list(X[:,k*blocksize:].todense().std(0).A1)
    X_std = np.array(X_std)
    return X_std


#   irlbpy code
#       From:   https://github.com/airysen/irlbpy
#       Date:   2021-12
#       License: Apache License, V 2.0, January 2004
#
#       Code unmodified except:
#           Added two print (feedback) blocks
#           Ran through lint formatting tool 'black' https://github.com/psf/black
#
import numpy as np
import scipy.sparse as sparse
import warnings

from numpy.fft import rfft, irfft
import numpy.linalg as nla


# Matrix-vector product wrapper
# A is a numpy 2d array or matrix, or a scipy matrix or sparse matrix.
# x is a numpy vector only.
# Compute A.dot(x) if t is False,  A.transpose().dot(x)  otherwise.


def multA(A, x, TP=False, L=None):
    if sparse.issparse(A):
        # m = A.shape[0]
        # n = A.shape[1]
        if TP:
            return sparse.csr_matrix(x).dot(A).transpose().todense().A[:, 0]
        return A.dot(sparse.csr_matrix(x).transpose()).todense().A[:, 0]
    if TP:
        return x.dot(A)
    return A.dot(x)


def multS(s, v, L, TP=False):
    N = s.shape[0]
    vp = prepare_v(v, N, L, TP=TP)
    p = irfft(rfft(vp) * rfft(s))
    if not TP:
        return p[:L]
    return p[L - 1 :]


def prepare_s(s, L=None):
    N = s.shape[0]
    if L is None:
        L = N // 2
    K = N - L + 1
    return np.roll(s, K - 1)


def prepare_v(v, N, L, TP=False):
    v = v.flatten()[::-1]
    K = N - L + 1
    if TP:
        lencheck = L
        if v.shape[0] != lencheck:
            raise VectorLengthException(
                "Length of v must be  L (if transpose flag is True)"
            )
        pw = K - 1
        v = np.pad(v, (pw, 0), mode="constant", constant_values=0)
    elif not TP:
        lencheck = N - L + 1
        if v.shape[0] != lencheck:
            raise VectorLengthException("Length of v must be N-K+1")
        pw = L - 1
        v = np.pad(v, (0, pw), mode="constant", constant_values=0)
    return v


def orthog(Y, X):
    """Orthogonalize a vector or matrix Y against the columns of the matrix X.
    This function requires that the column dimension of Y is less than X and
    that Y and X have the same number of rows.
    """
    dotY = multA(X, Y, TP=True)
    return Y - multA(X, dotY)


# Simple utility function used to check linear dependencies during computation:


def invcheck(x):
    # eps2 = 2 * np.finfo(np.float).eps
    eps2 = 2 * np.finfo(float).eps
    if x > eps2:
        x = 1 / x
    else:
        x = 0
        warnings.warn("Ill-conditioning encountered, result accuracy may be poor")
    return x


def lanczos(A, nval, tol=0.0001, maxit=50, center=None, scale=None, L=None):
    """Estimate a few of the largest singular values and corresponding singular
    vectors of matrix using the implicitly restarted Lanczos bidiagonalization
    method of Baglama and Reichel, see:

    Augmented Implicitly Restarted Lanczos Bidiagonalization Methods,
    J. Baglama and L. Reichel, SIAM J. Sci. Comput. 2005

    Keyword arguments:
    tol   -- An estimation tolerance. Smaller means more accurate estimates.
    maxit -- Maximum number of Lanczos iterations allowed.

    Given an input matrix A of dimension j * k, and an input desired number
    of singular values n, the function returns a tuple X with five entries:

    X[0] A j * nu matrix of estimated left singular vectors.
    X[1] A vector of length nu of estimated singular values.
    X[2] A k * nu matrix of estimated right singular vectors.
    X[3] The number of Lanczos iterations run.
    X[4] The number of matrix-vector products run.

    The algorithm estimates the truncated singular value decomposition:
    A.dot(X[2]) = X[0]*X[1].
    """

    import sys

    print(
        f">> lanczos A={A.shape}, nval={nval}, tol={tol}, maxit={maxit}",
        file=sys.stdout,
    )
    center_story = "None" if center is None else f"{center.shape}"
    scale_story = "None" if scale is None else f"{scale.shape}"
    print(f"++ lanczos center={center_story}, scale={scale_story}", file=sys.stdout)
    sys.stdout.flush()

    mmult = None
    m = None
    n = None
    if A.ndim == 2:
        mmult = multA
        m = A.shape[0]
        n = A.shape[1]
        if min(m, n) < 2:
            raise MatrixShapeException("The input matrix must be at least 2x2.")

    elif A.ndim == 1:
        mmult = multS
        A = np.pad(A, (0, A.shape[0] % 2), mode="edge")
        N = A.shape[0]
        if L is None:
            L = N // 2
        K = N - L + 1
        m = L
        n = K
        A = prepare_s(A, L)
    elif A.ndim > 2:
        raise MatrixShapeException("The input matrix must be 2D array")
    nu = nval

    m_b = min((nu + 20, 3 * nu, n))  # Working dimension size
    mprod = 0
    it = 0
    j = 0
    k = nu
    smax = 1
    # sparse = sparse.issparse(A)

    V = np.zeros((n, m_b))
    W = np.zeros((m, m_b))
    F = np.zeros((n, 1))
    B = np.zeros((m_b, m_b))

    V[:, 0] = np.random.randn(n)  # Initial vector
    V[:, 0] = V[:, 0] / np.linalg.norm(V)

    while it < maxit:
        if it > 0:
            j = k

        VJ = V[:, j]

        # apply scaling
        if scale is not None:
            VJ = VJ / scale

        W[:, j] = mmult(A, VJ, L=L)
        mprod = mprod + 1

        # apply centering
        # R code: W[, j_w] <- W[, j_w] - ds * drop(cross(dv, VJ)) * du
        if center is not None:
            W[:, j] = W[:, j] - np.dot(center, VJ)

        if it > 0:
            # NB W[:,0:j] selects columns 0,1,...,j-1
            W[:, j] = orthog(W[:, j], W[:, 0:j])
        s = np.linalg.norm(W[:, j])
        sinv = invcheck(s)
        W[:, j] = sinv * W[:, j]

        # Lanczos process
        while j < m_b:
            F = mmult(A, W[:, j], TP=True, L=L)
            mprod = mprod + 1

            # apply scaling
            if scale is not None:
                F = F / scale

            F = F - s * V[:, j]
            F = orthog(F, V[:, 0 : j + 1])
            fn = np.linalg.norm(F)
            fninv = invcheck(fn)
            F = fninv * F
            if j < m_b - 1:
                V[:, j + 1] = F
                B[j, j] = s
                B[j, j + 1] = fn
                VJp1 = V[:, j + 1]

                # apply scaling
                if scale is not None:
                    VJp1 = VJp1 / scale

                W[:, j + 1] = mmult(A, VJp1, L=L)
                mprod = mprod + 1

                # apply centering
                # R code: W[, jp1_w] <- W[, jp1_w] - ds * drop(cross(dv, VJP1))
                # * du
                if center is not None:
                    W[:, j + 1] = W[:, j + 1] - np.dot(center, VJp1)

                # One step of classical Gram-Schmidt...
                W[:, j + 1] = W[:, j + 1] - fn * W[:, j]
                # ...with full reorthogonalization
                W[:, j + 1] = orthog(W[:, j + 1], W[:, 0 : (j + 1)])
                s = np.linalg.norm(W[:, j + 1])
                sinv = invcheck(s)
                W[:, j + 1] = sinv * W[:, j + 1]
            else:
                B[j, j] = s
            j = j + 1
        # End of Lanczos process
        S = nla.svd(B)
        R = fn * S[0][m_b - 1, :]  # Residuals
        if it == 0:
            smax = S[1][0]  # Largest Ritz value
        else:
            smax = max((S[1][0], smax))

        conv = sum(np.abs(R[0:nu]) < tol * smax)
        if conv < nu:  # Not coverged yet
            k = max(conv + nu, k)
            k = min(k, m_b - 3)
        else:
            break
        # Update the Ritz vectors
        V[:, 0:k] = V[:, 0:m_b].dot(S[2].transpose()[:, 0:k])
        V[:, k] = F
        B = np.zeros((m_b, m_b))
        # Improve this! There must be better way to assign diagonal...
        for l in range(k):
            B[l, l] = S[1][l]
        B[0:k, k] = R[0:k]
        # Update the left approximate singular vectors
        W[:, 0:k] = W[:, 0:m_b].dot(S[0][:, 0:k])
        it = it + 1

    U = W[:, 0:m_b].dot(S[0][:, 0:nu])
    V = V[:, 0:m_b].dot(S[2].transpose()[:, 0:nu])
    # return((U, S[1][0:nu], V, it, mprod))

    print(f"<< lanczos it={it} mprod={mprod}", file=sys.stdout)
    sys.stdout.flush()

    return LanczosResult(
        **{"U": U, "s": S[1][0:nu], "V": V, "steps": it, "nmult": mprod}
    )


class LanczosResult:
    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])


class VectorLengthException(Exception):
    pass


class MatrixShapeException(Exception):
    pass


dataset_path = '/research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/datasets/'
analysis_save_path = '/research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/analysis/01-process/'

#adata = sc.read_h5ad(os.path.join(dataset_path, "Parse_10M_PBMC_cytokines.h5ad"))

#sc.pp.normalize_total(adata, target_sum=1e4)
#sc.pp.log1p(adata)
#sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.25)
#sc.pl.highly_variable_genes(adata) 
#adata.raw = adata
#adata = adata[:, adata.var.highly_variable].copy()

use_hv = True
#adata.uns['pca'] = {}
#adata.uns['pca']['params'] = {
#    'zero_center': True,
#    'use_highly_variable': use_hv,
#}

#adata.var['mean'] = adata.X.mean(0).A1
#adata.var['std'] = sparse_std(adata.X)

#S = lanczos(adata.X,50,center=adata.var['mean'],scale=adata.var['std'])

#adata.obsm['X_pca'] = (S.U * S.s)
#adata.varm['PCs'] = S.V
#adata.uns['pca']['variance'] = adata.obsm['X_pca'].var(0)
#adata.uns['pca']['variance_ratio'] = adata.uns['pca']['variance'] / (adata.X.shape[1] - 1 )

#adata.write(os.path.join(analysis_save_path, 'adata_post_pca.h5ad'))

adata = sc.read_h5ad(os.path.join(analysis_save_path, 'adata_post_pca.h5ad'))

sce.pp.harmony_integrate(adata, 'donor')
adata.write(os.path.join(analysis_save_path, 'adata_post_harmony.h5ad'))

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony")
sc.tl.umap(adata)
adata.write(os.path.join(analysis_save_path, 'adata_post_harmony_post_umap.h5ad'))

sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2)
adata.write(os.path.join(analysis_save_path, 'adata_post_harmony_post_umap_post_leiden.h5ad'))

UMAP_path = '/research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/analysis/01-process/UMAP/'
os.makedirs(UMAP_path, exist_ok=True)  # 如果目录不存在，则创建

# 绘制和保存 UMAP 图片
sc.pl.umap(adata, color='donor', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_donor.png"))
sc.pl.umap(adata, color='donor', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_donor.pdf"))

sc.pl.umap(adata, color='leiden', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_leiden.png"))
sc.pl.umap(adata, color='leiden', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_leiden.pdf"))

sc.pl.umap(adata, color='leiden', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_leiden_ondata.png"))
sc.pl.umap(adata, color='leiden', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_leiden_ondata.pdf"))

sc.pl.umap(adata, color='cell_type', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_cell_type.png"))
sc.pl.umap(adata, color='cell_type', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_cell_type.pdf"))

sc.pl.umap(adata, color='cell_type', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_cell_type_ondata.png"))
sc.pl.umap(adata, color='cell_type', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_cell_type_ondata.pdf"))

sc.pl.umap(adata, color='cytokine', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_cytokine.png"))
sc.pl.umap(adata, color='cytokine', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_cytokine.pdf"))

sc.pl.umap(adata, color='cytokine', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_cytokine_ondata.png"))
sc.pl.umap(adata, color='cytokine', legend_fontsize=8, legend_loc="on data", save=os.path.join(UMAP_path, "umap_cytokine_ondata.pdf"))

# 在 UMAP 图中显示 OLAH 基因的表达量
sc.pl.umap(adata, color='OLAH', cmap='viridis', legend_fontsize=8, save=os.path.join(UMAP_path, "umap_OLAH.png"))


### DEG
# 用于存储每个 cell_type 的差异基因信息
diff_gene_info = {}

# 遍历每种 cell_type
for cell_type in adata.obs['cell_type'].unique():
    # 子集化数据以包含当前 cell_type 的细胞
    adata_sub = adata[adata.obs['cell_type'] == cell_type]
    
    # 确保 cytokine 列为分类类型
    adata_sub.obs['cytokine'] = adata_sub.obs['cytokine'].astype("category")

    # 将 'PBS' 设置为参考组，并保持其他类别不变
    if 'PBS' in adata_sub.obs['cytokine'].cat.categories:
        adata_sub.obs['cytokine'].cat.set_categories(
            ['PBS'] + [cat for cat in adata_sub.obs['cytokine'].cat.categories if cat != 'PBS'],
            ordered=True,
            inplace=True
        )

    # 差异基因分析，PBS 作为参考
    sc.tl.rank_genes_groups(adata_sub, groupby='cytokine', reference='PBS', method='wilcoxon')

    # 提取差异基因信息并存储
    result = adata_sub.uns['rank_genes_groups']
    groups = result['names'].dtype.names  # 获取 cytokine 组名
    
    # 将所有信息存储到 DataFrame 中，并转换为便于读取的格式
    gene_data = {}
    for group in groups:
        gene_data[group] = pd.DataFrame({
            'gene': result['names'][group],
            'logfoldchanges': result['logfoldchanges'][group],
            'pvals': result['pvals'][group],
            'pvals_adj': result['pvals_adj'][group]
        })
    
    # 保存到主字典中
    diff_gene_info[cell_type] = gene_data


### Save and export

## Pickle format
import pickle
with open('diff_gene_info.pkl', 'wb') as f:
    pickle.dump(diff_gene_info, f)

## csv format
# 创建保存目录
csv_path = '/research/sharedresources/immunoinformatics/common/jqu/parse_biosciences/analysis/01-process/DEG_csv/'
os.makedirs(csv_path, exist_ok=True)

# 遍历字典，保存每个 DataFrame 为 CSV 文件
for cell_type, cytokine_dict in diff_gene_info.items():
    for cytokine, df in cytokine_dict.items():
        # 创建文件名，使用 {cell_type}_{cytokine}.csv 格式
        file_name = f"{cell_type}_{cytokine}.csv"
        # 将 DataFrame 保存为 CSV 文件
        df.to_csv(os.path.join(csv_path, file_name), index=True)

print("CSV files saved successfully.")


### Draw a heat map or matrix
# Column: cell_type
# Row: cytokine
# Cell: everage of OLAH

#import scanpy as sc
#import pandas as pd
#import numpy as np

# 假设 adata 是您的 AnnData 对象，并且已经运行了 Harmony 整合

# 1. 提取 OLAH 基因的表达数据（使用 Harmony 之后的值）
ol_expression = adata[:, 'OLAH'].X  # 获取 OLAH 基因的表达值
ol_expression_df = pd.DataFrame(ol_expression, index=adata.obs.index, columns=['OLAH'])

# 2. 将表达数据与 obs 中的 cell_type 和 cytokine 合并
ol_expression_df['cell_type'] = adata.obs['cell_type'].values
ol_expression_df['cytokine'] = adata.obs['cytokine'].values

with open('OLAH-1-expression.pkl', 'wb') as f:
    pickle.dump(ol_expression_df, f)

# 3. 计算每个 cell_type 和 cytokine 的平均表达量
heatmap_data = ol_expression_df.groupby(['cytokine', 'cell_type']).mean().unstack()

# 4. 用 NaN 填充缺失的值，确保热图中的所有单元格都有值
heatmap_data = heatmap_data.fillna(0)
with open('OLAH-2-heatmap_data.pkl', 'wb') as f:
    pickle.dump(heatmap_data, f)

# 5. 绘制热图
plt.figure(figsize=(10, 20))  # 设置图形大小
sns.heatmap(heatmap_data, cmap='viridis', annot=False)  # annot=True 显示每个单元格的值
plt.title('Mean Expression of OLAH (Post-Harmony)')
plt.xlabel('Cell Type')
plt.ylabel('Cytokine')

# 6. 保存热图为 PNG 和 PDF 格式
plt.savefig(os.path.join(analysis_save_path, 'mean_expression_OLAH.png'), bbox_inches='tight', dpi=300)
plt.savefig(os.path.join(analysis_save_path, 'mean_expression_OLAH.pdf'), bbox_inches='tight')  # 保存为 PDF


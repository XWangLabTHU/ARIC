"""
2021.5.21

@author: Wei Zhang, Hanwen Xu

E-mail: w-zhang16@mails.tsinghua.edu.cn, hw-xu16@mails.tsinghua.edu.cn

"""

from setuptools import setup, find_packages

setup(
    name="ARIC",
    version="0.0.1",
    keywords=("pip", "ARIC"),
    description=
    "ARIC: Accurate and robust inference of cell type proportions from bulk gene expression or DNA methylation data",
    license="GPL V3",
    url="",
    author="Wei Zhang, Hanwen Xu",
    author_email=
    "w-zhang16@mails.tsinghua.edu.cn, xuhw20@mails.tsinghua.edu.cn",
    packages=find_packages(),
    include_package_data=True,
    platforms="any",
    install_requires=[
        'numpy', 'tqdm', 'scipy', 'sklearn', 'statsmodels', 'pandas'
    ])

from setuptools import setup

setup(
    name='PyRINEX',
    version='3.0',
    author='Han Jinzhen',
    author_email='hanjinzhen9@gmail.com',

    long_description_content_type='text/markdown',
    url='https://github.com/yourusername/PyRINEX',
    packages=['PyRINEX'],
    classifiers=[
        "License :: OSI Approved :: Apache Software License 2.0",
    ],
    install_requires=[
        'numpy==1.21.1',
        'matplotlib==3.4.1',
        'chardet==4.0.0',
        'panda ==1.21.1',
        'seaborn==0.13.0'
    ],
    python_requires='>=3.6',
    license='Apache License 2.0'
)
wheel: deps
	NPY_DISABLE_CPU_FEATURES="AVX512F AVX512_SKX" python setup.py build bdist_wheel

deps:
	pip install wheel numpy setuptools

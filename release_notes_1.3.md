# BuildAMol meets Jax
BuildAMol can now use Google's Jax library to run some of its computations on GPU infrastructure. A new `use_jax()` global switch has been added, and some functions now support a `backend` argument which can be set to `jax` or `numba` (Numba support remains unchanged). A new `backends` package also facilitates integrating support for other libraries in the future such as PyTorch (not a supported backend at this time).


### Example
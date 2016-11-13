Run "demo.m"

To use your own masks, edit 'makeMeasurements.m'

The following files are key steps in the Polarization algorithm:

	makeMeasurements.m - makes noisy measurements using specified masks plus additional masks.
	reconstructPolarized.m - reconstructs signal via the polarization algorithm.
	relativeOptimalPhase.m - finds global shift constant.

The following files are used to demonstrate the Polarization algorithm:

	demo.m
	testPolarization.m

The following files are written by David F. Gleich (Stanford Univeristy) to find the largest component of a matrix:

	largest_component.m
	scomponents.m
	sparse_to_csr.m

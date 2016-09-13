import numpy as np
import astropy.table as tb
import astropy.io.fits as pf


def add_col(rec, name, arr=[], dtype=None):
	"""Generic function to add a new column to a structured numpy array.
	Borrows heavily from Tomek's code."""

	if len(arr)==0:
		arr=np.zeros(len(rec))

	arr = np.asarray(arr)
	if dtype is None:
		dtype = arr.dtype

	newdtype = np.dtype(rec.dtype.descr + [(name, dtype)])
	newrec = np.empty(rec.shape, dtype=newdtype)
	for field in rec.dtype.fields:
		newrec[field] = rec[field]

	newrec[name] = arr

	return newrec

def split_by(array, column_name, pivot, return_mask=False, logger=None):
	"""Return halves of a structured array split about a given value in one of the columns."""

	lower= array[(array[column_name]<pivot)]
	upper= array[(array[column_name]>pivot)]

	f = 1.0*len(upper)/len(array)

	if logger is not None:
		logger.info("Splitting array about %s=%3.4f (%3.2f percent above)"%(column_name,pivot, f*100.))

	if return_mask:
		return f, (array[column_name]<pivot)

	else:
		return f, upper, lower



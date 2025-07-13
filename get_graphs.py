import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d


def main():
	signal = pd.read_csv('signal.csv')
	signal.columns.tolist()
	dft = pd.read_csv('dft_result.csv')
	dft.columns.tolist()
	fft = pd.read_csv('fft_result.csv')
	fft.columns.tolist()

	try:
		time_dft = sys.argv[1]
	except:
		time_dft = 'N/A'
	
	try:
		time_fft = sys.argv[2]
	except:
		time_fft = 'N/A'

	fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8))

	ax1.plot(signal['x_axis'], signal['val'], 'b-o', label='Samples', markersize=3)
	'''x_new = np.linspace(signal['sample'].min(), signal['sample'].max(), 100)
	f_cubic = interp1d(signal['sample'], signal['val'], kind='cubic')
	y_linear = f_cubic(x_new)
	ax1.plot(x_new, y_linear, 'b--', label='Linear Interpolation', alpha=0.6)'''
	ax1.set_title('Sampled Signal')
	ax1.set_xlabel('sec [s]'),
	ax1.set_ylabel('f[k]')
	ax1.legend()
	ax1.grid(True)

	ax2.plot(dft['x_axis'], dft['val'], 'g-o', label='Samples', markersize=3)
	'''x_new = np.linspace(dft['sample'].min(), dft['sample'].max(), 100)
	f_cubic = interp1d(dft['sample'], dft['val'], kind='cubic')
	y_linear = f_cubic(x_new)
	ax2.plot(x_new, y_linear, 'g--', label='Linear Interpolation', alpha=0.6)'''
	ax2.set_title(f'DFT result - {time_dft} ms')
	ax2.set_xlabel('freq [Hz]')
	ax2.set_ylabel('mag(F[k])')
	ax2.legend()
	ax2.grid(True)

	ax3.plot(fft['x_axis'], fft['val'], 'r-o', label='Samples', markersize=3)
	'''x_new = np.linspace(fft['sample'].min(), fft['sample'].max(), 100)
	f_cubic = interp1d(fft['sample'], fft['val'], kind='cubic')
	y_linear = f_cubic(x_new)
	ax3.plot(x_new, y_linear, 'r--', label='Linear Interpolation', alpha=0.6)'''
	try:
		ax3.set_title(f'FFT result - {time_fft} ms - speed ratio: {float(time_dft)/float(time_fft):.3}')
	except:
		ax3.set_title(f'FFT result - {time_fft} ms')
	ax3.set_xlabel('freq [Hz]')
	ax3.set_ylabel('mag(F[k])')
	ax3.legend()
	ax3.grid(True)

	plt.tight_layout()
	plt.show()


if __name__ == '__main__':
	main()
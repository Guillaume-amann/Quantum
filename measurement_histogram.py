import pandas as pd
import matplotlib.pyplot as plt

hist = pd.read_csv("measurement_histogram.csv")
plt.plot(hist["state"], hist["count"])
plt.title("Measurement Histogram of Ground State")
plt.xlabel("Computational Basis State")
plt.ylabel("Counts (out of 1000)")
plt.tight_layout()
plt.show()
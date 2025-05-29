import matplotlib.pyplot as plt
import csv

# Replace this with your filename
FILENAME = "verify.csv"

s = []
value = []
fitted = []

with open(FILENAME, "r") as f:
    reader = csv.reader(f)
    header = next(reader)  # Skip header
    header, slope, intercept, error = (
        header[0:3],
        float(header[3]),
        float(header[4]),
        float(header[5]),
    )

    for row in reader:
        s.append(float(row[0]))
        value.append(float(row[1]))
        fitted.append(float(row[2]))

# Generate fitted line using slope and intercept
fitted_line = [slope * x + intercept for x in s]

# Plot
plt.figure(figsize=(8, 5))
plt.plot(s, value, "o", label="Original Value", color="blue")
plt.plot(s, fitted, "-", label="Fitted Value", color="green")
plt.plot(
    s, fitted_line, "--", label=f"Line: y = {slope:.2f}x + {intercept:.2f}", color="red"
)

plt.xlabel("s")
plt.ylabel("Value")
plt.title("Fitting Plot: " + str(error))
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

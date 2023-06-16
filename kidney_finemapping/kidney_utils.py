import matplotlib.pyplot as plt

# order of cell types in plot based on targets df
PLOT_KIDNEY_TARGET_INDICES = [2, 4, 8, 9, 7, 3, 6, 1, 5, 0]

# non-tubule epithelial cell types
NON_TUBULE_EPITHELIAL_TARGETS = ["Str", "End", "Tcell", "Immune"]

# colors for plotting
KIDNEY_CMAP = plt.get_cmap('tab10')
ROC_COLORS = {"PT": "#007d34",
              "Pod": "#817066",
              "CFH": "#803e75",
              "LOH": "#cea262",
              "DT": "#ff6800",
              "CD": "#ffb300"}

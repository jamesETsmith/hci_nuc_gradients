import os
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid


def make_figure_for_species(
    data_path: str, output_name: str, output_dir: str = "_figures"
):
    # active_space_list = [f"{N}e_{N}o" for N in [10, 20, 30, 40]]
    active_space_list = [f"{N}e_{N}o" for N in [40]]

    # Get smaller images to combine in a grid
    image_list = []
    for f in active_space_list:
        path = os.path.join(data_path, f)
        for i in ["118.png", "119.png", "120.png", "121.png"]:
            img_path = os.path.join(path, i)
            image_list.append(plt.imread(img_path))

    # Make figure
    fig = plt.figure()
    grid = ImageGrid(
        fig,
        111,  # similar to subplot(111)
        nrows_ncols=(1, 4),  # creates 3x4 grid of axes
        axes_pad=0.1,  # pad between axes in inch.
    )

    for ax, im in zip(grid, image_list):
        # Iterating over the grid returns the Axes.
        ax.imshow(im)
        ax.axis("off")

    # Save figure
    os.makedirs(output_dir, exist_ok=True)
    filename = os.path.join(output_dir, output_name)
    plt.tight_layout()
    plt.savefig(filename, dpi=900)


if __name__ == "__main__":

    make_figure_for_species("1_A", "HCISCF_NOs_1_A.png")
    make_figure_for_species("1_B", "HCISCF_NOs_1_B.png")
    make_figure_for_species("1_C", "HCISCF_NOs_1_C.png")

    make_figure_for_species("3_A", "HCISCF_NOs_3_A.png")
    make_figure_for_species("3_B", "HCISCF_NOs_3_B.png")
    make_figure_for_species("3_C", "HCISCF_NOs_3_C.png")

import open3d as o3d
import numpy as np
import os

def write_ply_from_pcd(points, colors, filename):
    vertex = np.zeros(points.shape[0], dtype=[("x", "f4"), ("y", "f4"), ("z", "f4"), ("red", "u1"), ("green", "u1"), ("blue", "u1")])
    vertex["x"] = points[:, 0]
    vertex["y"] = points[:, 1]
    vertex["z"] = points[:, 2]
    vertex["red"] = (colors[:, 0] * 255).astype(np.uint8)
    vertex["green"] = (colors[:, 1] * 255).astype(np.uint8)
    vertex["blue"] = (colors[:, 2] * 255).astype(np.uint8)

    with open(filename, 'w') as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("element vertex {}\n".format(len(points)))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property uchar red\n")
        f.write("property uchar green\n")
        f.write("property uchar blue\n")
        f.write("end_header\n")
        for i in range(len(points)):
            f.write("{:.3f} {:.3f} {:.3f} {} {} {}\n".format(points[i][0], points[i][1], points[i][2], vertex["red"][i], vertex["green"][i], vertex["blue"][i]))

def split_ply_by_color(merged_ply_path, output_directory):
    merged_ply = o3d.io.read_point_cloud(merged_ply_path)
    points = np.asarray(merged_ply.points)
    colors = np.asarray(merged_ply.colors)
    j = 0
    unique_colors, indices = np.unique(colors, axis=0, return_inverse=True)
    for i, color in enumerate(unique_colors):
        j = j +1
        color_points = points[indices == i]
        color_colors = colors[indices == i]
        color_str = "_".join(map(lambda x: str(int(x * 255)), color))  # Convert color to string for filename
        output_filename = os.path.join(output_directory, f"segment_{j}.ply")
        write_ply_from_pcd(color_points, color_colors, output_filename)

# Path to the merged PLY file
Merged_ply_path = r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_final/PLY_SDK_EU_Real_Merged_PCM_Translated_Rotated1.ply"

# Output directory for the separated PLY files
Output_directory = 'C:/Python Output/Real Data/PLY files/SDK_real_dataset_final/Segments/'
os.makedirs(Output_directory, exist_ok=True)

# Split the merged PLY file by color and save the separate PLY files
split_ply_by_color(Merged_ply_path, Output_directory)
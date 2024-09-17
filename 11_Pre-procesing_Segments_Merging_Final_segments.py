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

input_Segement_directory = 'C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/Segments/Final_Real_Segments/Unique_Colours/'
merged_wrong_point_cloud = o3d.geometry.PointCloud()

all_points = []
all_colors = []

for file_name_with_extension in os.listdir(input_Segement_directory):
    if file_name_with_extension.endswith('.ply'):
        filename, extension = os.path.splitext(file_name_with_extension)
        file_name_parts = filename.split('_')
        output_filename = f"{file_name_parts[0]}_{file_name_parts[1]}_{file_name_parts[2]}_{file_name_parts[3]}_{file_name_parts[4]}_{file_name_parts[5]}_{file_name_parts[6]}"
        filepath = os.path.join(input_Segement_directory, file_name_with_extension)
        Segment_ply = o3d.io.read_point_cloud(filepath)
        
        all_points.append(np.asarray(Segment_ply.points))
        all_colors.append(np.asarray(Segment_ply.colors))

all_points = np.concatenate(all_points, axis=0)
all_colors = np.concatenate(all_colors, axis=0)

Output_directory = os.path.join(input_Segement_directory, "Merged_Segment")
os.makedirs(Output_directory, exist_ok=True)
Merged_ply_path = os.path.join(Output_directory, output_filename + ".ply")

write_ply_from_pcd(all_points, all_colors, Merged_ply_path)
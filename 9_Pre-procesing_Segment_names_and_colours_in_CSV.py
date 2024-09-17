import open3d as o3d
import numpy as np
import os
import pandas as pd
def write_ply_from_pcd(points, colors, filename):
    vertex = np.zeros(points.shape[0], dtype=[("x", "f4"), ("y", "f4"), ("z", "f4"), ("red", "u1"), ("green", "u1"), ("blue", "u1")])
    vertex["x"] = points[:, 0]
    vertex["y"] = points[:, 1]
    vertex["z"] = points[:, 2]
    vertex["red"] = colors[:, 0]
    vertex["green"] = colors[:, 1]
    vertex["blue"] = colors[:, 2]

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
            f.write("{:.3f} {:.3f} {:.3f} {:.0f} {:.0f} {:.0f}\n".format(points[i][0], points[i][1], points[i][2], colors[i][0], colors[i][1], colors[i][2]))
Segm_Colour_name_array = np.array([["Segment Name", "Segment Colour" ]], dtype=object)
input_Segement_directory = 'C:/Python Output/Real Data/PLY files/SDK_real_dataset_1/Segments/Final_Real_Segments/Change their colour'
merged_wrong_point_cloud = o3d.geometry.PointCloud()
for file_name_with_extension in os.listdir(input_Segement_directory):
    if file_name_with_extension.endswith('.ply'):
        filename, extension = os.path.splitext(file_name_with_extension)
        
        filepath = os.path.join(input_Segement_directory, file_name_with_extension)
        Segment_ply = o3d.io.read_point_cloud(filepath)
        points = np.asarray(Segment_ply.points)
        colors = np.asarray(Segment_ply.colors)
        New_object = np.array([filename, colors[1]*255], dtype=object)
        Segm_Colour_name_array = np.vstack([Segm_Colour_name_array, New_object])

print(Segm_Colour_name_array)
def convert_to_csv(faces_array, csv_file_path):
    # Convert numpy array to pandas DataFrame
    df = pd.DataFrame(Segm_Colour_name_array[1:], columns=Segm_Colour_name_array[0])
    # Write to CSV
    df.to_csv(csv_file_path, index=False)

csv_file_pathXZ = r"C:/Python Output/Real data/CSV files/Names_Colours.csv"
convert_to_csv(Segm_Colour_name_array, csv_file_pathXZ)
#B.2.2 Process for super noisy point clouds
#this python script takes as input a super noisy ply file and aims to denoise it by creating 3D boxes 
import open3d as o3d
import numpy as np
import os 

pcd_file_path = r"C:/Python Output/Real data/PLY files/SDK_real_dataset_final/PLY_SDK_EU_Real.ply"

pcd = o3d.io.read_point_cloud(pcd_file_path)

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
# Convert point cloud to NumPy array for easier manipulation
points = np.asarray(pcd.points)
colors = np.asarray(pcd.colors)
print(points,colors)
# Ask user for U1 coordinates
coordinates_input1 = input("Enter x, y coordinates of upper point U1 separated by spaces: ")
# Split the input and convert to floats
x1, y1 = map(float, coordinates_input1.split())
# Create arrays to store the coordinates
X_coords= np.array([x1])
Y_coords= np.array([y1])
# Create an array U1 to store the coordinates
u1 = np.array([x1, y1])

# Ask user for U2 coordinates
coordinates_input2 = input("Enter x, y coordinates of upper point U2 separated by spaces: ")
# Split the input and convert to floats
x2, y2 = map(float, coordinates_input2.split())
# Add the new coordinates to arrays 
X_coords= np.hstack((X_coords, x2))
Y_coords= np.hstack((Y_coords, y2))
# Create an array U2 to store the coordinates
u2 = np.array([x2, y2])

# Ask user for U3 coordinates
coordinates_input3 = input("Enter x, y coordinates of upper point U3 separated by spaces: ")
# Split the input and convert to floats
x3, y3 = map(float, coordinates_input3.split())
# Add the new coordinates to arrays 
X_coords= np.hstack((X_coords, x3))
Y_coords= np.hstack((Y_coords, y3))
# Create an array U3 to store the coordinates
u3 = np.array([x3, y3])
# Ask user for U4 coordinates
coordinates_input4 = input("Enter x, y coordinates of upper point U4 separated by spaces: ")
# Split the input and convert to floats
x4, y4 = map(float, coordinates_input4.split())
# Add the new coordinates to arrays 
X_coords= np.hstack((X_coords, x4))
Y_coords= np.hstack((Y_coords, y4))
print("X coordinates:",X_coords,"Y coordinates:", Y_coords)
# Create an array U4 to store the coordinates
u4 = np.array([x4, y4])

z_max = float(input("Give cooridnate z of the upper boundary of the box Zmax: "))
z_minus = float(input("Give cooridnate z of the down boundary of the box Zmin: "))
# Find min and max values for each coordinate
x_minus = np.min(X_coords)
x_max = np.max(X_coords)
y_minus = np.min(Y_coords)
y_max = np.max(Y_coords)
# Select points within the cube boundaries
selected_indices = np.where(
    (points[:, 0] >= x_minus) & (points[:, 0] <= x_max) &
    (points[:, 1] >= y_minus) & (points[:, 1] <= y_max) &
    (points[:, 2] >= z_minus) & (points[:, 2] <= z_max)
)[0]

selected_points = points[selected_indices]
selected_colors = colors[selected_indices]

# Create a new point cloud with selected points and colors
selected_point_cloud = o3d.geometry.PointCloud()
selected_point_cloud.points = o3d.utility.Vector3dVector(selected_points)
selected_point_cloud.colors = o3d.utility.Vector3dVector(selected_colors)
# Split the file path into directory, filename, and extension
directory, file_name_with_extension = os.path.split(pcd_file_path)
file_name, extension = os.path.splitext(file_name_with_extension)
additional_word = "_Cube"
# Modify the filename based on the condition
new_file_name = file_name + additional_word
# Recreate the output file paths
pcd_file_path_out = [os.path.join(directory, new_file_name + extension)]
write_ply_from_pcd(selected_points, selected_colors*255, pcd_file_path_out[0] )
# Visualize original and translated point clouds
o3d.visualization.draw_geometries([selected_point_cloud])

